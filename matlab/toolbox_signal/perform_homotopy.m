function [X,Lambda] = perform_homotopy(D,y)

% perform_homotopy - compute the L1 regularization path
%
%   X = perform_homotopy(D,y);
%
%   Copyright (c) 2012 Gabriel Peyre

[P,N] = size(D);
niter = 10*P;

X = []; Lambda = [];
% initialization
C = D'*y;
[lambda,I] = max(abs(C));
x = zeros(N,1);
X(:,end+1) = x; Lambda(end+1) = lambda;
for it=1:niter
    % complementary support
    J = setdiff(1:N, I);
    % correlations
    c = D'*(y-D*x);
    if norm(c)<1e-6
        break;
    end
    % update direction
    d = zeros(N,1);
    d(I) = (D(:,I)'*D(:,I)) \ sign(c(I));
    % useful vector
    v = D(:,I)*d(I);
    % Compute minimum |gamma| so that situation 1) is in force.
    w = ( lambda-c(J) ) ./ ( 1 - D(:,J)'*v );
    gamma1 = min(w(w>0));
    if not(isempty(gamma1))
        i1 = J( w==gamma1 );
    end
    % Compute minimum |gamma| so that situation 2) is in force.
    w = ( lambda+c(J) ) ./ ( 1 + D(:,J)'*v );
    gamma2 = min(w(w>0));
    if not(isempty(gamma2))
        i2 = J( w==gamma2 );
    end
    % Compute minimum |gamma| so that situation 3) is in force.
    w = -x(I)./d(I);
    gamma3 = min(w(w>0));
    if not(isempty(gamma3))
        i3 = I( w==gamma3 );
    end
    % any condition is in force
    gamma = min([gamma1 gamma2 gamma3]);
    if isempty(gamma)
        break;  % final solution reached
    end
    % new solution
    x = x + gamma*d;
    lambda = lambda - gamma;
    if gamma==gamma1
        I = [I i1];
    elseif gamma==gamma2
        I = [I i2];
    elseif gamma==gamma3
        I(I==i3) = [];
        x(i3) = 0;
    end
    % record sparsity and lambda
    X(:,end+1) = x; Lambda(end+1) = lambda;
end


return;


n = size(D,2);

niter = getoptions(options, 'niter', min(round(size(D)))/2);

lambdaStop = 0;
solFreq = 1;
verbose = 0;
[X, numIters, activationHist] = SolveLasso(D, y, n, 'lasso', niter, lambdaStop, solFreq, verbose);
lambda_list = [];
sparsity_list = sum(X~=0);

return;



function [sols, numIters, activationHist] = SolveLasso(A, y, N, algType, maxIters, lambdaStop, solFreq, verbose, OptTol)
% SolveLasso: Implements the LARS/Lasso algorithms
% Usage
%	[sols, numIters, activationHist] = SolveLasso(A, y, N, algType,
%	maxIters, lambdaStop, solFreq, verbose, OptTol)
% Input
%	A           Either an explicit nxN matrix, with rank(A) = min(N,n) 
%               by assumption, or a string containing the name of a 
%               function implementing an implicit matrix (see below for 
%               details on the format of the function).
%	y           vector of length n.
%   N           length of solution vector. 
%   algType     'lars' for the Lars algorithm, 
%               'lasso' for lars with the lasso modification (default)
%	maxIters    maximum number of LARS iterations to perform. If not
%               specified, runs to stopping condition (default)
%   lambdaStop  If specified (and > 0), the algorithm terminates when the
%               Lagrange multiplier <= lambdaStop. 
%   solFreq     if =0 returns only the final solution, if >0, returns an 
%               array of solutions, one every solFreq iterations (default 0). 
%   verbose     1 to print out detailed progress at each iteration, 0 for
%               no output (default)
%	OptTol      Error tolerance, default 1e-5
% Outputs
%	sols           solution of the Lasso/LARS problem
%	numIters       Total number of steps taken
%   activationHist Array of indices showing elements entering and 
%                  leaving the solution set
% Description
%   SolveLasso implements the LARS algorithm, as described by Efron et al. in 
%   "Least Angle Regression". Currently, the original algorithm is
%   implemented, as well as the lasso modification, which solves 
%      min || y - Ax ||_2^2 s.t || x ||_1 <= t
%   for all t > 0, using a path following method with parameter t.
%   The implementation implicitly factors the active set matrix A(:,I)
%   using Choleskly updates. 
%   The matrix A can be either an explicit matrix, or an implicit operator
%   implemented as an m-file. If using the implicit form, the user should
%   provide the name of a function of the following format:
%     y = OperatorName(mode, m, n, x, I, dim)
%   This function gets as input a vector x and an index set I, and returns
%   y = A(:,I)*x if mode = 1, or y = A(:,I)'*x if mode = 2. 
%   A is the m by dim implicit matrix implemented by the function. I is a
%   subset of the columns of A, i.e. a subset of 1:dim of length n. x is a
%   vector of length n is mode = 1, or a vector of length m is mode = 2.
% References
%   B. Efron, T. Hastie, I. Johnstone and R. Tibshirani, 
%   "Least Angle Regression", Annals of Statistics, 32, 407-499, 2004
% See Also
%   SolveOMP, SolveBP, SolveStOMP
%

if nargin < 9,
    OptTol = 1e-5;
end
if nargin < 8,
    verbose = 0;
end
if nargin < 7,
    solFreq = 0;
end
if nargin < 6,
    lambdaStop = 0;
end
if nargin < 5,
    maxIters = 10*N;
end
if nargin < 4,
    algType = 'lasso';
end

switch lower(algType)
    case 'lars'
        isLasso = 0;
    case 'lasso'
        isLasso = 1;
end

explicitA = ~(ischar(A) || isa(A, 'function_handle'));
n = length(y);

% Global variables for linsolve function
global opts opts_tr machPrec
opts.UT = true; 
opts_tr.UT = true; opts_tr.TRANSA = true;
machPrec = 1e-5;

x = zeros(N,1);
iter = 0;

% First vector to enter the active set is the one with maximum correlation
if (explicitA)
    corr = A'*y;             
else
    corr = feval(A,2,n,N,y,1:N,N); % = A'*y             
end
lambda = max(abs(corr));
newIndices = find(abs(abs(corr)-lambda) < machPrec)';    

% Initialize Cholesky factor of A_I
R_I = [];
activeSet = [];
for j = 1:length(newIndices)
    iter = iter+1;
    [R_I, flag] = updateChol(R_I, n, N, A, explicitA, activeSet, newIndices(j));
    activeSet = [activeSet newIndices(j)];
    if verbose
        fprintf('Iteration %d: Adding variable %d\n', iter, activeSet(j));
    end
end

activationHist = activeSet;
collinearIndices = [];
sols = [];

done = 0;
while  ~done
    % Compute Lars direction - Equiangular vector
    dx = zeros(N,1);
    % Solve the equation (A_I'*A_I)dx_I = corr_I
    z = linsolve(R_I,corr(activeSet),opts_tr);
    dx(activeSet) = linsolve(R_I,z,opts);
    if (explicitA)
        dmu = A'*(A(:,activeSet)*dx(activeSet));
    else
        dmu = feval(A,1,n,length(activeSet),dx(activeSet),activeSet,N); 
        dmu = feval(A,2,n,N,dmu,1:N,N); 
    end

    % For Lasso, Find first active vector to violate sign constraint
    if isLasso
        zc = -x(activeSet)./dx(activeSet);
        gammaI = min([zc(zc > 0); inf]);
        removeIndices = activeSet(find(zc == gammaI));
    else
        gammaI = Inf;
        removeIndices = [];
    end

    % Find first inactive vector to enter the active set
    if (length(activeSet) >= min(n, N))
        gammaIc = 1;
    else
        inactiveSet = 1:N;
        inactiveSet(activeSet) = 0;
        inactiveSet(collinearIndices) = 0;
        inactiveSet = find(inactiveSet > 0);
        lambda = abs(corr(activeSet(1)));
        dmu_Ic = dmu(inactiveSet);
        c_Ic = corr(inactiveSet);

        epsilon = 1e-12; 
        gammaArr = [(lambda-c_Ic)./(lambda - dmu_Ic + epsilon) (lambda+c_Ic)./(lambda + dmu_Ic + epsilon)]';
        gammaArr(gammaArr < machPrec) = inf;
        gammaArr = min(gammaArr)';
        [gammaIc, Imin] = min(gammaArr);
    end

    % If gammaIc = 1, we are at the LS solution
    if (1-gammaIc) < OptTol
        newIndices = [];
    else
        newIndices = inactiveSet(find(abs(gammaArr - gammaIc) < machPrec));
        %newIndices = inactiveSet(Imin);
    end

    gammaMin = min(gammaIc,gammaI);

    % Compute the next Lars step
    x = x + gammaMin*dx;
    corr = corr - gammaMin*dmu;

    % Check stopping condition
    if ((1-gammaMin) < OptTol) | ((lambdaStop > 0) & (lambda <= lambdaStop))
        done = 1;
    end

    % Add new indices to active set
    if (gammaIc <= gammaI) && (length(newIndices) > 0)
        for j = 1:length(newIndices)
            iter = iter+1;
            if verbose
                fprintf('Iteration %d: Adding variable %d\n', iter, newIndices(j));
            end
            % Update the Cholesky factorization of A_I
            [R_I, flag] = updateChol(R_I, n, N, A, explicitA, activeSet, newIndices(j));
            % Check for collinearity
            if (flag)
                collinearIndices = [collinearIndices newIndices(j)];
                if verbose
                    fprintf('Iteration %d: Variable %d is collinear\n', iter, newIndices(j));
                end
            else
                activeSet = [activeSet newIndices(j)];
                activationHist = [activationHist newIndices(j)];
            end
        end
    end

    % Remove violating indices from active set
    if (gammaI <= gammaIc)
        for j = 1:length(removeIndices)
            iter = iter+1;
            col = find(activeSet == removeIndices(j));
            if verbose
                fprintf('Iteration %d: Dropping variable %d\n', iter, removeIndices(j));
            end
            % Downdate the Cholesky factorization of A_I
            R_I = downdateChol(R_I,col);
            activeSet = [activeSet(1:col-1), activeSet(col+1:length(activeSet))];
            
            % Reset collinear set
            collinearIndices = [];
        end

        x(removeIndices) = 0;  % To avoid numerical errors
        activationHist = [activationHist -removeIndices];
    end

    if iter >= maxIters
        done = 1;
    end

    if done | ((solFreq > 0) & (~mod(iter,solFreq)))
        sols = [sols x];
    end
end

numIters = iter;
clear opts opts_tr machPrec


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R, flag] = updateChol(R, n, N, A, explicitA, activeSet, newIndex)
% updateChol: Updates the Cholesky factor R of the matrix 
% A(:,activeSet)'*A(:,activeSet) by adding A(:,newIndex)
% If the candidate column is in the span of the existing 
% active set, R is not updated, and flag is set to 1.

global opts_tr machPrec
flag = 0;

if (explicitA)
    newVec = A(:,newIndex);
else
    e = zeros(N,1);
    e(newIndex) = 1;
    newVec = feval(A,1,n,N,e,1:N,N); 
end

if length(activeSet) == 0,
    R = sqrt(sum(newVec.^2));
else
    if (explicitA)
        p = linsolve(R,A(:,activeSet)'*A(:,newIndex),opts_tr);
    else
        AnewVec = feval(A,2,n,length(activeSet),newVec,activeSet,N);
        p = linsolve(R,AnewVec,opts_tr);
    end
    q = sum(newVec.^2) - sum(p.^2);
    if (q <= machPrec) % Collinear vector
        flag = 1;
    else
        R = [R p; zeros(1, size(R,2)) sqrt(q)];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = downdateChol(R, j)
% downdateChol: `Downdates' the cholesky factor R by removing the 
% column indexed by j.

% Remove the j-th column
R(:,j) = [];
[m,n] = size(R);

% R now has nonzeros below the diagonal in columns j through n.
% We use plane rotations to zero the 'violating nonzeros'.
for k = j:n
    p = k:k+1;
    [G,R(p,k)] = planerot(R(p,k));
    if k < n
        R(p,k+1:n) = G*R(p,k+1:n);
    end
end

% Remove last row of zeros from R
R = R(1:n,:);


%
% Copyright (c) 2006. Yaakov Tsaig and Joshua Sweetkind-Singer
%  

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
