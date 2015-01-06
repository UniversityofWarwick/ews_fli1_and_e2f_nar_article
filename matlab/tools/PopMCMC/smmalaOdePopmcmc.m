function [ xPosterior, BurnInPath, logLikelihood ] = smmalaOdePopmcmc(model, varargin)

%========== Parsing Sampler Parameters ===========
samplerP = inputParser;
samplerP.addParamValue('Iterations', 1000, @(x) x > 0 );
samplerP.addParamValue('BurnIn', 500, @(x) x >= 0 );
samplerP.addParamValue('TemperatureDensity', @(x) (x./max(x)).^4, @(x) isa(x,'function_handle') );
samplerP.addParamValue('StepSize', 0.1);
samplerP.addParamValue('Verbose', 1000);
samplerP.addParamValue('Adapt', 100);
samplerP.addParamValue('LogFile', 1);
samplerP.addParamValue('ThinBurnIn',1);
samplerP.addParamValue('Init','model');
samplerP.parse( varargin{:} );

NumOfIterations = samplerP.Results.Iterations;
BurnIn =  samplerP.Results.BurnIn;
TemperatureDensity = samplerP.Results.TemperatureDensity;
StepSize = samplerP.Results.StepSize;
Verbose = samplerP.Results.Verbose;
Adapt = samplerP.Results.Adapt;
LogFile = samplerP.Results.LogFile;
ThinBurnIn = samplerP.Results.ThinBurnIn;
Init = samplerP.Results.Init;

clear samplerP;


if ischar(LogFile)
    if numlabs > 1
        [path, name, ext] = fileparts(LogFile);
        LogFile = fullfile(path, [name '_' num2str(labindex) ext]);
    end
    LogFile = fopen(LogFile, 'w');
end
%=================================================


%Message IDs and broadcast node
MPI_MSG_LL = 10;
MPI_MSG_STATE = 20;
MPI_BROADCAST_ID = 1;

%======= Pre-allocations and initilisation =======
D = model.NumberOfParameters;
xPosterior = zeros(NumOfIterations - BurnIn, D);
if nargout > 2
    logLikelihood = zeros(NumOfIterations - BurnIn,1);
end

if nargout > 1
    bunInLength = floor((BurnIn-1)/ThinBurnIn) + 1;
    BurnInPath = zeros(bunInLength, D);
end

if ischar(Init)
    X = model.SamplePrior(1);
else
    X = Init;
end

tp = 1:numlabs;
t = TemperatureDensity(tp);

[LL dLL ML] = model.LLCmp(X);
[LP dLP MP] = model.PLCmp(X);

StepSize_par = StepSize*t(labindex);
propPar = 0;
accPar = 0;

PropLad = 0;
AccLad = 0;
if labindex == MPI_BROADCAST_ID
    propEx = zeros(numlabs,1);
    accEx  = zeros(numlabs,1);
end
%=================================================

for IterationNum = 1:NumOfIterations
    %============== Main Sampling Loop ===============
    
    %=== Update Power Posteriors ====
    % This is done in paralel.
    e = StepSize_par;
    propPar = propPar + 1;
    
    x_b  = X;
    
    %====== Simplified MMALA proposal ========
    dL_b = dLL*t(labindex) + dLP;
    G_b  = t(labindex)^2*ML + MP;
    cholG_b = chol(G_b);
    
    mean_pb = x_b + (e/2)*(cholG_b\(cholG_b'\dL_b));
    new_x = mean_pb + 1/sqrt(e)*cholG_b\randn(D,1);
    
    pno = sum(log(diag(cholG_b))) - 0.5*sum( (cholG_b*(new_x - mean_pb)).^2 )/e;
    
    %calculate probability of old state from the new state
    [new_LL new_dLL new_ML] = model.LLCmp(new_x);
    [new_LP new_dLP new_MP] = model.PLCmp(new_x);
    
    dL_b = new_dLL*t(labindex) + new_dLP;
    G_b  = t(labindex)^2*new_ML + new_MP;
    [cholG_b ec] = chol(G_b);
    
    if ~ec
        % if G is positive definite continue otherwie reject.
        
        mean_pb = new_x + (e/2)*(cholG_b\(cholG_b'\dL_b));
        
        pon = sum(log(diag(cholG_b))) - 0.5*sum( (cholG_b*(x_b - mean_pb)).^2 )/e;
        
        %Accept/reject
        JL_b  = t(labindex)*LL + LP;
        new_JL_b = t(labindex)*new_LL + new_LP;
        
        Ratio = new_JL_b + pon - JL_b - pno;
        
        if Ratio > 0 || (Ratio > log(rand))
            accPar = accPar +1;
            X = new_x;
            LL = new_LL;
            LP = new_LP;
            dLL = new_dLL;
            dLP = new_dLP;
            ML = new_ML;
            MP = new_MP;
        end
        
    end
    
    % =========== Random exchanges between chains ===========
    
    if numlabs > 1
        if labindex == MPI_BROADCAST_ID
            for i = 1:numlabs
                PropLad = PropLad + 1;
                % select a random pair of nodes i,j
                idx_i = randi(numlabs);
                switch idx_i
                    case 1
                        idx_j = 2;
                    case numlabs
                        idx_j = numlabs-1;
                    otherwise
                        if rand > 0.5
                            idx_j = idx_i+1;
                        else
                            idx_j = idx_i-1;
                        end
                end
                
                propEx(idx_i) = propEx(idx_i) + 1;
                propEx(idx_j) = propEx(idx_j) + 1;
                
                %fprintf('Pair of nodes (%g,%g) \n',idx_i,idx_j);
                
                % Get LL from nodes i and j
                %fprintf('Request LL from nodes %g %g \n',idx_i,idx_j);
                labBroadcast(MPI_BROADCAST_ID, {idx_i, idx_j, 'LL'});
                if idx_i ~= MPI_BROADCAST_ID && idx_j ~= MPI_BROADCAST_ID
                    LL_i = labReceive(idx_i, MPI_MSG_LL);
                    %fprintf('LL received from node %g \n',idx_i);
                    LL_j = labReceive(idx_j, MPI_MSG_LL);
                    %fprintf('LL received from node %g \n',idx_j);
                else
                    if idx_i == MPI_BROADCAST_ID
                        LL_i = LL;
                        LL_j = labReceive(idx_j, MPI_MSG_LL);
                        %fprintf('LL received from node %g \n',idx_j);
                    else
                        LL_i = labReceive(idx_i, MPI_MSG_LL);
                        %fprintf('LL received from node %g \n',idx_i);
                        LL_j = LL;
                    end
                end
                
                % Calculate acceptance ratio for exchange move
                Ratio = LL_j*t(idx_i) + LL_i*t(idx_j) - LL_i*t(idx_i) - LL_j*t(idx_j);
                
                if Ratio > 0 || (Ratio > log(rand))
                    AccLad = AccLad + 1;
                    accEx(idx_i) = accEx(idx_i) + 1;
                    accEx(idx_j) = accEx(idx_j) + 1;
                    % Notify nodes i and j to exchange states.
                    %fprintf('Notify node %g to exchange with node %g\n',idx_i,idx_j)
                    labBroadcast(MPI_BROADCAST_ID, {idx_i, idx_j, 'state'});
                    
                    if idx_i == MPI_BROADCAST_ID
                        %fprintf('node %g waits to exchange with node %g\n',idx_i,idx_j)
                        State_j = labSendReceive(idx_j, idx_j, X, MPI_MSG_STATE);
                        %fprintf('node %g exchanged with node %g\n',idx_i,idx_j)
                        X = State_j;
                        [LL dLL ML] = model.LLCmp(X);
                        [LP dLP MP] = model.PLCmp(X);
                    end
                    if idx_j == MPI_BROADCAST_ID
                        %fprintf('node %g waits to exchange with node %g\n',idx_j, idx_i)
                        State_i = labSendReceive(idx_i, idx_i, X, MPI_MSG_STATE);
                        %fprintf('node %g exchanged with node %g\n',idx_j, idx_i);
                        X = State_i;
                        [LL dLL ML] = model.LLCmp(X);
                        [LP dLP MP] = model.PLCmp(X);
                    end
                end
            end
            % Notify nodes to exit sync loop.
            %fprintf('Notify nodes to exit sync loop\n');
            labBroadcast(MPI_BROADCAST_ID, 'exit');
            
        else
            % all nodes except MPI_BROADCAST_ID enter into sync loop.
            while (true)
                %Get broadcasted messages and act
                message = labBroadcast(MPI_BROADCAST_ID);
                
                if ischar(message) && strcmpi(message,'exit')
                    %fprintf('Exit message received\n');
                    break; % exit sync loop.
                end
                if message{1} == labindex || message{2} == labindex
                    switch message{3}
                        case 'LL'
                            % send LL to MPI_BROADCAST_ID
                            %fprintf('Node %g sending LL \n', labindex);
                            labSend(LL, MPI_BROADCAST_ID, MPI_MSG_LL);
                        case 'state'
                            % exchange state with toLab
                            if message{1} == labindex
                                toLab = message{2};
                            else
                                toLab = message{1};
                            end
                            %fprintf('Node %g waits exchange with node %g\n',labindex,toLab)
                            State_new = labSendReceive(toLab,toLab, X, MPI_MSG_STATE);
                            %fprintf('Node %g exchanged with node %g\n',labindex,toLab)
                            X = State_new;
                            [LL dLL ML] = model.LLCmp(X);
                            [LP dLP MP] = model.PLCmp(X);
                            
                        otherwise
                            fprintf(LogFile,'node %g received unknown message \n',labindex);
                    end
                end
            end
        end
    end
    
    % ========= Sampler Bookkeeping and Logging  ===========
    
    if mod(IterationNum, Verbose) == 0
        if labindex == MPI_BROADCAST_ID
            fprintf(LogFile,'Itteration: %g \n', IterationNum);
            fprintf(LogFile,'Exchange rate: %g \n',AccLad/PropLad);
            for chain_i = 1:numlabs
                fprintf(LogFile,'\t == Exchange rate for chain %g: %g \n', chain_i, accEx(chain_i)/propEx(chain_i));
            end
        end
        fprintf(LogFile,'Node: %g, Temp.: %g, Acc.: %g \n', labindex, t(labindex), accPar/propPar );
    end
    
    %Step-Size tunning during Burnin phase
    if (mod(IterationNum, Adapt) == 0) && (IterationNum < BurnIn)
        %fprintf(1,'Tunning Stepsizes ...\n');
        
        AccRatio = accPar/propPar;
        
        if AccRatio > 0.75
            StepSize_par = StepSize_par*1.2;
        elseif AccRatio < 0.68
            StepSize_par = StepSize_par*0.8;
        end
        fprintf(LogFile,'Node: %g, Temp.: %g, StepSize: %g \n',labindex, t(labindex),StepSize_par);
        accPar = 0;
        propPar = 0;
    end
    
    
    % Save samples
    if IterationNum > BurnIn
        xPosterior(IterationNum-BurnIn,:) = X;
        if nargout > 2
            logLikelihood(IterationNum-BurnIn) = LL;
        end
    end
    
    if nargout > 1 && IterationNum <= BurnIn && mod( (IterationNum-1), ThinBurnIn) == 0
        itThinIdx = floor((IterationNum-1)/ThinBurnIn) + 1;
        BurnInPath(itThinIdx,:) = X;
    end
    
    % Start timer after burn-in
    if IterationNum == BurnIn
        if labindex == MPI_BROADCAST_ID
            fprintf(LogFile,'Burn-in complete, now drawing posterior samples.\n');
        end
    end
    
end

end
