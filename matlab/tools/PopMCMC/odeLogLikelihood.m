function [ LL, dLL, FI ] = odeLogLikelihood( Ti, X, epsilon, theta, odeRhs, y0, tStart, A, thetaS, stateIdx)
% Ti: time points
% X: observations
% odeRhs: compiled library
% A: transformation matrix
% thetaS: parameter selector for calculating steady state conditions
t0 = 0;
if nargin > 6 && ~isempty(tStart)
    t0 = tStart;
end

Ti    = Ti(:);
y0    = y0(:);
theta = theta(:);

No = size(X,1);
Np = size(theta,1);
Nt = length(Ti);

Ny = size(y0,1);

dY0 = [zeros(Ny,Np) eye(Ny)];
Np = Np+Ny;

%preallocate tmp matrices.
if nargin > 7 && ~isempty(A)
    kron_eyeNpA = sparse(kron(eye(Np),A));
end
if nargout > 1 
    eyeNpNo = sparse(kron(eye(Np),ones(1,No)));
    oneNp1 = ones(Np,1);
end

if nargin > 8 && ~isempty(thetaS)
    %%% initial conditions are subject to steady state solution.
    dY0 = [zeros(Ny,length(theta)) eye(Ny)];
    theta0 = thetaS.*theta;
    
    %empty second param solves until steady state.
    [yI, yIp] = odeSolverMex(odeRhs, [], theta0, y0, dY0, [], t0);
    
    yIp = kron([thetaS;ones(Ny,1)],ones(Ny,1)).*yIp;
    
    % replace initial conditions with the new estimates.
    y0 = yI;
    dY0 = reshape(yIp,[Ny Np]);
end

nTi = Ti;
if Ti(1) == t0     
    nTi = Ti(2:end);
end

if nargout > 1
    [y, yp] = odeSolverMex(odeRhs, nTi, theta, y0, dY0, [], t0);
else
    y = odeSolverMex(odeRhs, nTi, theta, y0, dY0, [], t0);
end

if Ti(1) == t0
    y = [y0 y];
    if nargout > 1
        yp = [dY0(:) yp];
    end
end

if nargin > 7 && ~isempty(A)
    y = A*y;
    if nargout > 1
        yp = kron_eyeNpA*yp;
    end
end

I = ~isnan(X);
X(isnan(X)) = 0;
I1 = I.*(X-y);
I2 = I1.*(X-y);

LL = -0.5*sum(sum(I2(stateIdx{1},:)))/epsilon(1) -0.5*sum(sum(I2(stateIdx{2},:)))/epsilon(2) ... 
    -0.5*log(2*pi*epsilon(1))*sum(sum(I(stateIdx{1},:))) -0.5*log(2*pi*epsilon(2))*sum(sum(I(stateIdx{2},:)));
% LL = -Nt*No*0.5*log(epsilon) - 0.5*sum(sum((X - y).^2,1))./epsilon;

% Gradient
if nargout > 1

    dLL = zeros(Np + 2,1);
    sI1 = I1;
    sI1(stateIdx{1},:) = sI1(stateIdx{1},:)/epsilon(1);
    sI1(stateIdx{2},:) = sI1(stateIdx{2},:)/epsilon(2);
    
    dLL(1:Np,1) = sum( eyeNpNo * (kron(oneNp1,sI1).*yp), 2);
    % dLL(1:Np,1) = sum( eyeNpNo * (kron(oneNp1,(X - y)./epsilon).*yp), 2);
   
    dLL(end-1,1)  = 0.5*sum(sum(I2(stateIdx{1},:)))/(epsilon(1)^2) -0.5*sum(sum(I(stateIdx{1},:)))/epsilon(1);
    dLL(end,1)  = 0.5*sum(sum(I2(stateIdx{2},:)))/(epsilon(2)^2) -0.5*sum(sum(I(stateIdx{2},:)))/epsilon(2);
    % dLL(end,1)  = -Nt*No*0.5./epsilon + 0.5*sum(sum((X - y).^2,1))./epsilon.^2;
    
    % Fisher information
    if nargout > 2
        
        Iyp = kron(oneNp1, I).*yp;
        FI = zeros(Np+2);
        
        Iyp_f =reshape(Iyp,[No Np Nt]);
        
        yPt_r = reshape( permute(Iyp_f(stateIdx{1},:), [2 1 3]), Np, length(stateIdx{1})*Nt);
        yPt_p = reshape( permute(Iyp_f(stateIdx{2},:), [2 1 3]), Np, length(stateIdx{2})*Nt);
        % yPt = reshape( permute( reshape(yp,[No Np Nt]), [2 1 3]), Np, No*Nt);
        FI(1:Np,1:Np) = (yPt_r*yPt_r')/epsilon(1) + (yPt_p*yPt_p')/epsilon(2);
        FI(end-1,end-1)   = 0.5*sum(sum(I(stateIdx{1},:)))/(epsilon(1)^2);
        FI(end-1,end-1)   = 0.5*sum(sum(I(stateIdx{2},:)))/(epsilon(2)^2);
        % FI(1:Np,1:Np) = (yPt*yPt')./epsilon ;
        % FI(end,end)   = Nt*No*0.5/epsilon.^2;
        
%         if numel(epsilon) == 1
%             FI = (yPt*yPt')./epsilon ;
%         else
%             yPtE = bsxfun(@(x,y) x./y, yPt, epsilon(:)');
%             FI = yPtE*yPt';
%         end
        
    end
end

end
