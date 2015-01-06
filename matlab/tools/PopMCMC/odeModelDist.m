classdef odeModelDist < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        
        %ode model related
        name;           % name for the odeModel
        Ny;             % number of ode system state variables
        Np;             % number of ode system parameters
        odeRHS;         % ode system compiled shared library
        parameters;     % array with the description of parameters
        states;         % array with the description of state variables
        
        %prior related
        priorStd;
        priorMean;
        
        mData;          
        tStart;         % Starting time for the integrator, T0
        mThetaS;        % Parameter selector for steady state initial 
                        % conditions.
    end
    
    properties (Dependent)
        Name;
        Data;
        Parameters;
        States;
        NumberOfParameters;
        NormalPrior;
    end
    
    methods (Static)
        
        function model = odeModelDistFromData(Data, odeSystem, solib, prior)
            % model = odeModelDistFromData(Data, odeSystem, prior)
            %   Creates an ode model object instance.
            %   Data - a structure or stucture array with fields:
            %   Data.Observations - 
            %        NxT matrix with measurments  of N variables 
            %        at T time points.
            %   Data.TimePoints   - 
            %        1xT row vector of T time points.
            %   Data.TransformationMatrix - 
            %        A NxV matrix mapping the V system variables 
            %        to the N observations. If all system variables are
            %        observed then this can be an empty array.
            %   odeSystem - a string with the vf file of the ode model.
            %   prior - a 1x2 cell array wiht prior(1) the mean for the
            %           of the prior and prior(2) the std.
            %           Both the mean and std can be vectors with size
            %           equal to the number of parameters.
            %
            %tmpLogFile02 = fopen('/home/ucaktpa/tmp/tmpLogFile02.txt', 'w');
            %fprintf(tmpLogFile02, 'Started evaluation of odeModelDistFromData\n')
            switch nargin
                case 3
                    %fprintf(tmpLogFile02, 'Started evaluation of case 2 of switch nargin in odeModelDistFromData\n')
                    %if exist(odeSystem, 'file')
                    %  fprintf(tmpLogFile02, 'Model VF file exists inside odeModelDistFromData\n')
                    %end
                    %if exist(odeSoLib, 'file')
                    %  fprintf(tmpLogFile02, 'Model dynamic librafy file exists inside odeModelDistFromData\n')
                    %end
                    %fclose(tmpLogFile02)
                    model = odeModelDist(odeSystem,solib);
                case 4 
                    model = odeModelDist(odeSystem,solib, prior);
            end
            
            model.Data = Data;
        end
        
    end
    
    methods

        function self = odeModelDist(odeSystemFile,solib, prior)
            %tmpLogFile03 = fopen('/home/ucaktpa/tmp/tmpLogFile03.txt', 'w');           
            %fprintf(tmpLogFile03, 'Started evaluation of odeModelDist(odeSystemFile)\n')
            %fclose(tmpLogFile03)
            %tmpLogFile04 = fopen('/home/ucaktpa/tmp/tmpLogFile04.txt', 'w');
            %fprintf(tmpLogFile04, 'Started evaluation of readXMLmodel(odeSystemFile)\n')
            self.readXMLmodel(odeSystemFile);
            %fprintf(tmpLogFile04, 'Completed evaluation of readXMLmodel(odeSystemFile)\n')
            %fclose(tmpLogFile04)
            %tmpLogFile05 = fopen('/home/ucaktpa/tmp/tmpLogFile05.txt', 'w');
            %fprintf(tmpLogFile05, 'Started assignment of odeSoLibFile\n')
            %self.odeRHS = compileModel(odeSystemFile);
            %system_vf =  which(odeSystemFile);

            %self.odeRHS =  compileModel(system_vf);
            self.odeRHS  = solib;
            %fprintf(tmpLogFile05, 'Completed assignment of odeSoLibFile\n')
            %fclose(tmpLogFile05)            

            self.tStart = 0;
            
            self.priorStd = 1;
            self.priorMean = 0;
            
            if nargin > 2 && ~isempty(prior)
                self.NormalPrior = prior;
            else
                self.priorStd = 1;
                self.priorMean = 0;
            end
            
        end
        
        function sim = simulateODE(self, timepoints, params)
            timepoints = timepoints(:);
            if nargin > 2 && ~isempty(params)
                x = params;
            else
                x = self.getDefaultParameters();
            end
            theta = 10.^x(1:self.Np);
            y0 = 10.^x(self.Np+1 : end-2);
            
            if ~isempty(self.mThetaS)
                theta0 = theta.*self.mThetaS;
                y0 = odeSolverMex(self.odeRHS, [], theta0, y0, [], [], self.tStart);
            end
            
            tps = timepoints;
            if timepoints(1) == self.tStart
                tps = timepoints(2:end,1);
            end
            
            sim = odeSolverMex(self.odeRHS, tps, theta, y0, [], [], self.tStart);
            
            if timepoints(1) == self.tStart
                sim = [y0 sim];
            end
            
        end
        
        function sim = simulateData(self, timepoints, params, transformationMatrix)
            if nargin > 2 && ~isempty(params)
                x = params;
            else
                x = self.getDefaultParameters();
            end
            
            sim = simulateODE(self, timepoints, x);
            
            if nargin > 3 && ~isempty(transformationMatrix);
                sim = transformationMatrix*simulateODE(self, timepoints, params);
            end
            sim = sim + normrnd(0,10.^x(end),size(sim));
        end
        
        function [LL dLL ML] = LLCmp(self, params)
            % LL logLikelihood, dLL gradient and ML Fisher information
            if isempty(self.mData)
                error('odeModelDist:NotInitialised','Add data using the Data property of this object');
            end
            
            x = 10.^params;
            theta = x(1 : self.Np);
            y0 = x(self.Np+1 : end-2);
            epsilon = x(end-1:end);
            
            LL  = 0;
            dLL = zeros(self.NumberOfParameters,1);
            FI  = zeros(self.NumberOfParameters);
            
            for i = 1:length(self.mData)
                [cLL, cdLL, cFI]=odeLogLikelihood(self.mData(i).TimePoints, ...
                                                  self.mData(i).Observations,...
                                                  epsilon, theta, self.odeRHS,...
                                                  y0, self.tStart,...
                                                  self.mData(i).TransformationMatrix,...
                                                  self.mThetaS, self.mData(i).StateIdx);
                LL = LL + cLL;
                dLL = dLL + cdLL;
                FI = FI + cFI;
            end
            
            if nargout > 1
                J = diag(x.*log(10));
                dLL = J*dLL;
                if nargout > 2
                    ML = J'*FI*J;
                end
            end
        end
        
        function [PL dPL MP] = PLCmp(self, params)
            % Prior term
            D = self.NumberOfParameters;
            % calculate prior
            PL = -0.5*sum(((params - self.priorMean).^2)./self.priorStd);
            
            if nargout > 1
                dPL = -(params - self.priorMean)./self.priorStd;
                if nargout > 2
                    MP  = diag(ones(1,D)./self.priorStd);
                end
            end
        end
        
        function [JLL dJLL MJL] = JLLCmp(self, params)
            % Joint log-likelihood = LLCmp+PLCmp
            switch nargout
                case 1
                    JLL = self.LLCmp(params) + self.PLCmp(params);
                case 2
                    [LL, dLL] = self.LLCmp(params);
                    [LP, dLP] = self.PLCmp(params);
                    JLL = LL + LP;
                    dJLL = dLL + dLP;
                case 3
                    [LL, dLL, FI] = self.LLCmp(params);
                    [LP, dLP, M] = self.PLCmp(params);
                    JLL = LL + LP;
                    dJLL = dLL + dLP;
                    MJL = FI + M;
                otherwise
                    error('odeModelDist:InvalidParameters','Too many parameters');
            end
        end
            
        function newX = SamplePrior(self, numOfSamples)
            % draw random sample from the prior
            if isscalar(self.priorStd)
                newX = normrnd(self.priorMean, self.priorStd, self.NumberOfParameters, numOfSamples);
            else
                newX = normrnd(self.priorMean, repmat(self.priorStd, 1, numOfSamples), numOfSamples);
            end
        end
        
        function ret = get.Data(self)
            ret = self.mData;
        end
        
        function set.Data(self, dataArray)
            
            for i = 1:length(dataArray)
                [e,msg] = self.isValidDataStructure(dataArray(i));
                if e 
                    error('odeModelDist:InvalidDataStructure',['Data(' num2str(i) '): ' msg]);
                end
            end
            self.mData = dataArray;
            %[self.mTimePoints, self.mTP2DataIDX] = self.getCommonTimePoints(dataArray);
            
        end
        
        function NumberOfParameters = get.NumberOfParameters(self)
            % unknown parameters: kinetic rates+ initial conditions + std
            %NumberOfParameters = self.Ny + self.Np + 1;
            NumberOfParameters = length(self.parameters);
        end
        
        function [MeanStd] = get.NormalPrior(self)
            MeanStd = cell(2,1);
            MeanStd{1} = self.priorMean;
            MeanStd{2} = self.priorStd;
        end
        
        function set.NormalPrior(self, MeanStd)
            if numel(MeanStd) ~= 2
                error('odeModelDist:NormalPrior:InvalidParaneters','prior should have two elements mean and std');
            end
            
            if iscell(MeanStd)
                if ~isscalar(MeanStd{1}) && numel(MeanStd{1}) ~= self.NumberOfParameters
                    error('odeModelDist:NormalPrior:InvalidParaneters','The mean should be a scalar or of the same dimensions as the number of parameters');
                end
                
                if ~isscalar(MeanStd{2}) && numel(MeanStd{2}) ~= self.NumberOfParameters
                    error('odeModelDist:NormalPrior:InvalidParaneters','Std should be a scalar or of the same dimension as the number of parameters');
                end
                
                self.priorMean = MeanStd{1};
                self.priorStd  = MeanStd{2};
            else
                self.priorMean = MeanStd(1);
                self.priorStd  = MeanStd(2);
            end
        end
        
        function [strName] = get.Name(self)
            strName = self.name;
        end
        
        function [ret] = get.States(self)
            ret = self.states;
        end
        
        function [ret] = get.Parameters(self)
            ret = self.parameters;
        end
        
        function [ret] = getDefaultParameters(self)
            ret = log10([self.parameters(:).DefaultValue]');
        end
    end
    
    methods (Access = private)
        
%         function [TP idx] = getCommonTimePoints(self, Data)
%             lengths = arrayfun(@(x) length(x.TimePoints), Data);    
%             idx = cell2mat(arrayfun(@(i,j) i*ones(1,j), 1:length(Data), lengths,'UniformOutput',false));
%             TP = cat(2, Data(:).TimePoints);
%             
%             [TP, isort] = sort(TP);
%             idx = idx(isort);
%             
%             iend = find(diff([TP Inf]));
%             istart = [1 iend(1:end-1)+1];
%             idx = arrayfun(@(i,j) idx(i:j), istart,iend,'UniformOutput',false);
%             TP = unique(TP);
%         end
            
        function readXMLmodel(self, xmlFile)
            %read xmlFile and set the states and parameters arrays,
            tree = xmlread(xmlFile);
            if ~tree.hasChildNodes
                error('odeModelDist:InvalidFileFormat','Invalid ode model file format: %s.',xmlFile);
            end
            
            cNodes = tree.getChildNodes;
            if cNodes.getLength ~= 1
                error('odeModelDist:InvalidFileFormat','Invalid ode model file format: %s.',xmlFile);
            end
            
            vfNode = cNodes.item(0);
            strName = char(vfNode.getAttribute('Name'));
            if isempty(strName)
                error('odeModelDist:InvalidFileFormat','Invalid ode model file format: %s.',xmlFile);
            end
            self.name = strName;
            
            stateNodes = vfNode.getElementsByTagName('StateVariable');
            numStateNodes = stateNodes.getLength;
            if numStateNodes < 1
                error('odeModelDist:InvalidFileFormat','Invalid ode model file format: %s. No state variables found.',xmlFile);
            end
            
            %self.states
            self.Ny = numStateNodes;
            for c = 1:numStateNodes
                stateNode = stateNodes.item(c-1);
                self.states(c).Name = char( stateNode.getAttribute('Name') );
                self.states(c).Description = char( stateNode.getAttribute('Description') );
                self.states(c).Formula =  char( stateNode.getAttribute('Formula') );
                self.states(c).DefaultInitialCondition =  str2double( stateNode.getAttribute('DefaultInitialCondition') );
            end
            
            
            paramNodes = vfNode.getElementsByTagName('Parameter');
            numParamNodes = paramNodes.getLength;
            if numParamNodes < 1
                error('odeModelDist:InvalidFileFormat','Invalid ode model file format: %s. No parameter variables found.',xmlFile);
            end
            self.Np = numParamNodes;
            
            % unknown parameters: kinetic rates+ initial conditions + std 
            for c = 1:numParamNodes
                paramNode = paramNodes.item(c-1);
                self.parameters(c).Name = char( paramNode.getAttribute('Name') );
                self.parameters(c).Description = char( paramNode.getAttribute('Description') );
                self.parameters(c).DefaultValue =  str2double( paramNode.getAttribute('DefaultValue') );
            end
            
            for c = 1:numStateNodes
                self.parameters(numParamNodes+c).Name = [self.states(c).Name '_0'];
                self.parameters(numParamNodes+c).Description = ['Initial conditions for state variable ' self.states(c).Name];
                self.parameters(numParamNodes+c).DefaultValue = self.states(c).DefaultInitialCondition;
            end
            
            self.parameters(numParamNodes+numStateNodes+1).Name = 'epsilon_r';
            self.parameters(numParamNodes+numStateNodes+1).Description = 'Standard deviation of observation noise for mRNA states';
            self.parameters(numParamNodes+numStateNodes+1).DefaultValue = 1;
            
            self.parameters(numParamNodes+numStateNodes+2).Name = 'epsilon_p';
            self.parameters(numParamNodes+numStateNodes+2).Description = 'Standard deviation of observation noise for Protein states';
            self.parameters(numParamNodes+numStateNodes+2).DefaultValue = 1;
            
            %Steady state parameters
            self.mThetaS = [];
            selectParamsNodes = vfNode.getElementsByTagName('SteadyStateParams');
            if selectParamsNodes.getLength > 0
                paramExludeStr = char ( selectParamsNodes.item(0).getAttribute('ExcludeList') );
                if isempty(paramExludeStr)
                    error('odeModelDist:InvalidFileFormat','Invalid ode model file format: %s. SteadyStateParams defined but empty.',xmlFile);
                end
                paramsExcludeLst = textscan(paramExludeStr,'%s','Delimiter',',');
                paramsExcludeLst  = paramsExcludeLst{1};
                if length(paramsExcludeLst)<1
                    error('odeModelDist:InvalidFileFormat','Invalid ode model file format: %s. SteadyStateParams invalid list format.',xmlFile);
                end
                
                self.mThetaS = ones(self.Np,1);
                
                paramNames =  {self.parameters(1:self.Np).Name};
                for i = 1:length(paramsExcludeLst);
                    idx = find(ismember(paramNames,paramsExcludeLst{i}),1);
                    if isempty(idx)
                        error('odeModelDist:InvalidFileFormat','Parameter %s specified in the SteadyStateParams section of file %s is not defined as a parameter.',paramsExcludeLst{i}, xmlFile);
                    end
                    self.mThetaS(idx) = 0;
                end
            end
            
        end
        
        
        function [err msg] = isValidDataStructure(self, Data)
            
            if ~isfield(Data,'Observations') || isempty(Data.Observations)
                msg = 'Could not find field Observations in Data';
                err = true;
                return;
            end
            
            if ~isfield(Data,'TimePoints') || isempty(Data.TimePoints)
                msg = 'Could not find field TimePoints in Data';
                err = true;
                return;
            end
            
            if size(Data.TimePoints,1) ~= 1
                msg = 'TimePoints must be a row vector';
                err = true;
                return
            end
            
            Ntps = size(Data.TimePoints,2);
            if any(diff(Data.TimePoints) <= 0)
                msg = 'TimePoints must be a monotonically increasing series of values';
                err = true;
                return
            end
            
            if size(Data.Observations,2) ~= Ntps
                msg = 'The second dimension of observations should match the size of timePoints';
                err = true;
                return;
            end
            
            Nobs = size(Data.Observations,1);
            if (self.Ny ~= Nobs) && (~isfield(Data,'TransformationMatrix') || isempty(Data.TransformationMatrix))
                msg = 'Size of ode system does not match the observations. You need to provide a tranformation matrix for the observation model.';
                err = true;
                return;
            end
            
            if isfield(Data,'TransformationMatrix') && ~isempty(Data.TransformationMatrix)
                [Ar, Ac] = size(Data.TransformationMatrix);
                if Ar ~= Nobs || Ac ~= self.Ny
                   msg = 'Size of transformationMatrix does not much the size of observation or system states';
                   err = true;
                   return;
                end
            end
            
            msg = '';
            err = false;
        end
        
    end
    
end
