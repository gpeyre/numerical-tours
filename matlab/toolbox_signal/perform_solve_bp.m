function sol = p(A, y, N, maxIters, lambda, OptTol)
% SolveBP: Solves a Basis Pursuit problem
% Usage
%	sol = SolveBP(A, y, N, maxIters, lambda, OptTol)
% Input
%	A           Either an explicit nxN matrix, with rank(A) = min(N,n)
%               by assumption, or a string containing the name of a
%               function implementing an implicit matrix (see below for
%               details on the format of the function).
%	y           vector of length n.
%   N           length of solution vector.
%	maxIters    maximum number of PDCO iterations to perform, default 20.
%   lambda      If 0 or omitted, Basis Pursuit is applied to the data,
%               otherwise, Basis Pursuit Denoising is applied with
%               parameter lambda (default 0).
%	OptTol      Error tolerance, default 1e-3
% Outputs
%	 sol        solution of BP
% Description
%   SolveBP solves the basis pursuit problem
%      min ||x||_1 s.t. A*x = y
%   by reducing it to a linear program, and calling PDCO, a primal-dual
%   log-barrier algorithm. Alternatively, if lambda ~= 0, it solves the
%   Basis Pursuit Denoising (BPDN) problem
%      min lambda*||x||_1 + 1/2||y - A*x||_2^2
%   by transforming it to an SOCP, and calling PDCO.
%   The matrix A can be either an explicit matrix, or an implicit operator
%   implemented as a function. If using the implicit form, the user should
%   provide the name of a function of the following format:
%     y = OperatorName(mode, m, n, x, I, dim)
%   This function gets as input a vector x and an index set I, and returns
%   y = A(:,I)*x if mode = 1, or y = A(:,I)'*x if mode = 2.
%   A is the m by dim implicit matrix implemented by the function. I is a
%   subset of the columns of A, i.e. a subset of 1:dim of length n. x is a
%   vector of length n is mode = 1, or a vector of length m is mode = 2.
% See Also
%   SolveLasso, SolveOMP, SolveITSP
%

if nargin < 6,
    OptTol = 1e-3;
end
if nargin < 5,
    lambda = 0;
end
if nargin < 4,
    maxIters = 20;
end

n = length(y);

n_pdco = 2*N;    % Input size
m_pdco = n;      % Output size

% upper and lower bounds
bl = zeros(n_pdco,1);
bu = Inf .* ones(n_pdco,1);

% generate the vector c
if (lambda ~= 0)
    c = lambda .* ones(n_pdco,1);
else
    c = ones(n_pdco,1);
end

% Generate an initial guess
x0 = ones(n_pdco,1)/n_pdco;       % Initial x
y0 = zeros(m_pdco,1);             % Initial y
z0 = ones(n_pdco,1)/n_pdco;       % Initial z

d1 = 1e-4;                 % Regularization parameters
if (lambda ~= 0) % BPDN
    d2 = 1;
else
    d2 = 1e-4;
end

xsize = 1;                 % Estimate of norm(x,inf) at solution
zsize = 1;                 % Estimate of norm(z,inf) at solution

options = pdcoSet;         % Option set for the function pdco
options = pdcoSet( options, ...
    'MaxIter    ', maxIters  , ...
    'FeaTol     ', OptTol    , ...
    'OptTol     ', OptTol    , ...
    'StepTol    ', 0.99      , ...
    'StepSame   ', 0         , ...
    'x0min      ', 0.1       , ...
    'z0min      ', 1.0       , ...
    'mu0        ', 0.01      , ...
    'method     ', 1         , ...
    'LSQRMaxIter', 20        , ...
    'LSQRatol1  ', 1e-3      , ...
    'LSQRatol2  ', 1e-15     , ...
    'wait       ', 0    );

if (ischar(A) || isa(A, 'function_handle'))
    [xx,yy,zz,inform,PDitns,CGitns,time] = ...
        pdco(c, @pdcoMat, y, bl, bu, d1, d2, options, x0, y0, z0, xsize, zsize);
else
    Phi = [A -A];
    [xx,yy,zz,inform,PDitns,CGitns,time] = ...
        pdco(c, Phi, y, bl, bu, d1, d2, options, x0, y0, z0, xsize, zsize);
end

% Extract the solution from the output vector x
sol = xx(1:N) - xx((N+1):(2*N));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = pdcoMat(mode,m,n,x)
if (mode == 1) % Direct operator
    % Decompose input
    n2 = n/2;
    u = x(1:n2);
    v = x(n2+1:n);
    
    % Apply matrix A
    Au = feval(A,1,m,n2,u,1:n2,n2);
    Av = feval(A,1,m,n2,v,1:n2,n2);
    
    y = Au-Av;
else % Adjoint operator
    n2 = n/2;
    Atx = feval(A,2,m,n2,x,1:n2,n2);
    y = [Atx; -Atx];
end


%
% Copyright (c) 2006. Yaakov Tsaig
%

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%

function [x,y,z,inform,PDitns,CGitns,time] = ...
    pdco( Fname,Aname,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize )

%-----------------------------------------------------------------------
% pdco.m: Primal-Dual Barrier Method for Convex Objectives (23 Sep 2003)
%-----------------------------------------------------------------------
% [x,y,z,inform,PDitns,CGitns,time] = ...
%   pdco(Fname,Aname,b,bl,bu,d1,d2,options,x0,y0,z0,xsize,zsize);
%
% solves optimization problems of the form
%
%    minimize    phi(x) + 1/2 norm(D1*x)^2 + 1/2 norm(r)^2
%      x,r
%    subject to  A*x + D2*r = b,   bl <= x <= bu,   r unconstrained,
%
% where
%    phi(x) is a smooth convex function  defined by function Fname;
%    A      is an m x n matrix defined by matrix or function Aname;
%    b      is a given m-vector;
%    D1, D2 are positive-definite diagonal matrices defined from d1, d2.
%           In particular, d2 indicates the accuracy required for
%           satisfying each row of Ax = b.
%
% D1 and D2 (via d1 and d2) provide primal and dual regularization
% respectively.  They ensure that the primal and dual solutions
% (x,r) and (y,z) are unique and bounded.
%
% A scalar d1 is equivalent to d1 = ones(n,1), D1 = diag(d1).
% A scalar d2 is equivalent to d2 = ones(m,1), D2 = diag(d2).
% Typically, d1 = d2 = 1e-4.
% These values perturb phi(x) only slightly  (by about 1e-8) and request
% that A*x = b be satisfied quite accurately (to about 1e-8).
% Set d1 = 1e-4, d2 = 1 for least-squares problems with bound constraints.
% The problem is then equivalent to
%
%    minimize    phi(x) + 1/2 norm(d1*x)^2 + 1/2 norm(A*x - b)^2
%    subject to  bl <= x <= bu.
%
% More generally, d1 and d2 may be n and m vectors containing any positive
% values (preferably not too small, and typically no larger than 1).
% Bigger elements of d1 and d2 improve the stability of the solver.
%
% At an optimal solution, if x(j) is on its lower or upper bound,
% the corresponding z(j) is positive or negative respectively.
% If x(j) is between its bounds, z(j) = 0.
% If bl(j) = bu(j), x(j) is fixed at that value and z(j) may have
% either sign.
%
% Also, r and y satisfy r = D2 y, so that Ax + D2^2 y = b.
% Thus if d2(i) = 1e-4, the i-th row of Ax = b will be satisfied to
% approximately 1e-8.  This determines how large d2(i) can safely be.
%
%
% EXTERNAL FUNCTIONS:
% options         = pdcoSet;                  provided with pdco.m
% [obj,grad,hess] = Fname( x );               provided by user
%               y = Aname( name,mode,m,n,x ); provided by user if pdMat
%                                             is a string, not a matrix
%
% INPUT ARGUMENTS:
% Fname      may be an explicit n x 1 column vector c,
%            or a string containing the name of a function Fname.m
%%%!!! Revised 12/16/04 !!!
%            (Fname cannot be a function handle)
%            such that [obj,grad,hess] = Fname(x) defines
%            obj  = phi(x)              : a scalar,
%            grad = gradient of phi(x)  : an n-vector,
%            hess = diag(Hessian of phi): an n-vector.
%         Examples:
%            If phi(x) is the linear function c'*x, Fname could be
%            be the vector c, or the name or handle of a function
%            that returns
%               [obj,grad,hess] = [c'*x, c, zeros(n,1)].
%            If phi(x) is the entropy function E(x) = sum x(j) log x(j),
%            Fname should return
%               [obj,grad,hess] = [E(x), log(x)+1, 1./x].
% Aname      may be an explicit m x n matrix A (preferably sparse!),
%            or a string containing the name of a function Aname.m
%%%!!! Revised 12/16/04 !!!
%            (Aname cannot be a function handle)
%            such that y = aname( name,mode,m,n,x )
%            returns   y = A*x (mode=1)  or  y = A'*x (mode=2).
%            The input parameter "name" will be the string 'Aname'
%            or whatever the name of the actual function is.
% b          is an m-vector.
% bl         is an n-vector of lower bounds.  Non-existent bounds
%            may be represented by bl(j) = -Inf or bl(j) <= -1e+20.
% bu         is an n-vector of upper bounds.  Non-existent bounds
%            may be represented by bu(j) =  Inf or bu(j) >=  1e+20.
% d1, d2     may be positive scalars or positive vectors (see above).
% options    is a structure that may be set and altered by pdcoSet
%            (type help pdcoSet).
% x0, y0, z0 provide an initial solution.
% xsize, zsize are estimates of the biggest x and z at the solution.
%            They are used to scale (x,y,z).  Good estimates
%            should improve the performance of the barrier method.
%
%
% OUTPUT ARGUMENTS:
% x          is the primal solution.
% y          is the dual solution associated with Ax + D2 r = b.
% z          is the dual solution associated with bl <= x <= bu.
% inform = 0 if a solution is found;
%        = 1 if too many iterations were required;
%        = 2 if the linesearch failed too often.
%        = 3 if the step lengths became too small.
% PDitns     is the number of Primal-Dual Barrier iterations required.
% CGitns     is the number of Conjugate-Gradient  iterations required
%            if an iterative solver is used (LSQR).
% time       is the cpu time used.
%----------------------------------------------------------------------

% PRIVATE FUNCTIONS:
%    pdxxxbounds
%    pdxxxdistrib
%    pdxxxlsqr
%    pdxxxlsqrmat
%    pdxxxmat
%    pdxxxmerit
%    pdxxxresid1
%    pdxxxresid2
%    pdxxxstep
%
% GLOBAL VARIABLES:
%    global pdDDD1 pdDDD2 pdDDD3
%
%
% NOTES:
% The matrix A should be reasonably well scaled: norm(A,inf) =~ 1.
% The vector b and objective phi(x) may be of any size, but ensure that
% xsize and zsize are reasonably close to norm(x,inf) and norm(z,inf)
% at the solution.
%
% The files defining Fname  and Aname
% must not be called Fname.m or Aname.m!!
%
%
% AUTHOR:
%    Michael Saunders, Systems Optimization Laboratory (SOL),
%    Stanford University, Stanford, California, USA.
%    saunders@stanford.edu
%
% CONTRIBUTORS:
%    Byunggyoo Kim, Samsung, Seoul, Korea.
%    hightree@samsung.co.kr
%
% DEVELOPMENT:
% 20 Jun 1997: Original version of pdsco.m derived from pdlp0.m.
% 29 Sep 2002: Original version of pdco.m  derived from pdsco.m.
%              Introduced D1, D2 in place of gamma*I, delta*I
%              and allowed for general bounds bl <= x <= bu.
% 06 Oct 2002: Allowed for fixed variabes: bl(j) = bu(j) for any j.
% 15 Oct 2002: Eliminated some work vectors (since m, n might be LARGE).
%              Modularized residuals, linesearch
% 16 Oct 2002: pdxxx..., pdDDD... names rationalized.
%              pdAAA eliminated (global copy of A).
%              Aname is now used directly as an explicit A or a function.
%              NOTE: If Aname is a function, it now has an extra parameter.
% 23 Oct 2002: Fname and Aname can now be function handles.
% 01 Nov 2002: Bug fixed in feval in pdxxxmat.
% 19 Apr 2003: Bug fixed in pdxxxbounds.
% 07 Aug 2003: Let d1, d2 be scalars if input that way.
% 10 Aug 2003: z isn't needed except at the end for output.
% 10 Aug 2003: mu0 is now an absolute value -- the initial mu.
% 13 Aug 2003: Access only z1(low) and z2(upp) everywhere.
%              stepxL, stepxU introduced to keep x within bounds.
%              (With poor starting points, dx may take x outside,
%              where phi(x) may not be defined.
%              Entropy once gave complex values for the gradient!)
% 16 Sep 2003: Fname can now be a vector c, implying a linear obj c'*x.
% 19 Sep 2003: Large system K4 dv = rhs implemented.
% 23 Sep 2003: Options LSproblem and LSmethod replaced by Method.
% 18 Nov 2003: stepxL, stepxU gave trouble on lptest (see 13 Aug 2003).
%              Disabled them for now.  Nonlinear problems need good x0.
% 19 Nov 2003: Bugs with x(fix) and z(fix).
%              In particular, x(fix) = bl(fix) throughout, so Objective
%              in iteration log is correct for LPs with explicit c vector.
%-----------------------------------------------------------------------

global pdDDD1 pdDDD2 pdDDD3

if 0
    myprintf('\n   --------------------------------------------------------')
    myprintf('\n   pdco.m                            Version of 19 Nov 2003')
    myprintf('\n   Primal-dual barrier method to minimize a convex function')
    myprintf('\n   subject to linear constraints Ax + r = b,  bl <= x <= bu')
    myprintf('\n   --------------------------------------------------------\n')
end

m = length(b);
n = length(bl);

%---------------------------------------------------------------------
% Decode Fname.
%---------------------------------------------------------------------
%%%!!! Revised 12/16/04 !!!
% Fname cannot be a function handle

operator  =  ischar(Fname);
explicitF = ~operator;

if explicitF
    myprintf('\n')
    mydisp('The objective is linear')
else
    fname = Fname;
    myprintf('\n')
    mydisp(['The objective function  is named       ' fname])
end

%---------------------------------------------------------------------
% Decode Aname.
%---------------------------------------------------------------------
%%%!!! Revised 12/16/04 !!!
% Aname cannot be a function handle
operator  =  ischar(Aname) || isa(Aname, 'function_handle');
explicitA = ~operator;

if explicitA            % assume Aname is an explicit matrix A.
    nnzA   = nnz(Aname);
    if  issparse(Aname)
        myprintf('The matrix A is an explicit sparse matrix')
    else
        myprintf('The matrix A is an explicit dense matrix' )
    end
    myprintf('\n\nm        = %8g     n        = %8g      nnz(A)  =%9g', m,n,nnzA)
else
    if ischar(Aname)
        mydisp(['The matrix A is an operator defined by ' Aname])
    end
    myprintf('\nm        = %8g     n        = %8g', m,n)
end

normb  = norm(b ,inf);   normx0 = norm(x0,inf);
normy0 = norm(y0,inf);   normz0 = norm(z0,inf);
myprintf('\nmax |b | = %8g     max |x0| = %8.1e', normb , normx0)
myprintf(                '      xsize   = %8.1e', xsize)
myprintf('\nmax |y0| = %8g     max |z0| = %8.1e', normy0, normz0)
myprintf(                '      zsize   = %8.1e', zsize)

%---------------------------------------------------------------------
% Initialize.
%---------------------------------------------------------------------
true   = 1;
false  = 0;
zn     = zeros(n,1);
nb     = n + m;
nkkt   = nb;
CGitns = 0;
inform = 0;
% 07 Aug 2003: No need for next lines.
%if length(d1)==1, d1 = d1*ones(n,1); end   % Allow scalar d1, d2
%if length(d2)==1, d2 = d2*ones(m,1); end   % to mean d1*e, d2*e

%---------------------------------------------------------------------
% Grab input options.
%---------------------------------------------------------------------
maxitn    = options.MaxIter;
featol    = options.FeaTol;
opttol    = options.OptTol;
steptol   = options.StepTol;
stepSame  = options.StepSame;   % 1 means stepx==stepz
x0min     = options.x0min;
z0min     = options.z0min;
mu0       = options.mu0;
Method    = options.Method;
itnlim    = options.LSQRMaxIter * min(m,n);
atol1     = options.LSQRatol1;  % Initial  atol
atol2     = options.LSQRatol2;  % Smallest atol, unless atol1 is smaller
conlim    = options.LSQRconlim;
wait      = options.wait;

%---------------------------------------------------------------------
% Set other parameters.
%---------------------------------------------------------------------
kminor    = 0;      % 1 stops after each iteration
eta       = 1e-4;   % Linesearch tolerance for "sufficient descent"
maxf      = 10;     % Linesearch backtrack limit (function evaluations)
maxfail   = 1;      % Linesearch failure limit (consecutive iterations)
bigcenter = 1e+3;   % mu is reduced if center < bigcenter
thresh    = 1e-8;   % For sparse LU with Method=41


% Parameters for LSQR.
atolmin   = eps;    % Smallest atol if linesearch back-tracks
btol      = 0;      % Should be small (zero is ok)
show      = false;  % Controls LSQR iteration log
gamma     = max(d1);
delta     = max(d2);

myprintf('\n\nx0min    = %8g     featol   = %8.1e', x0min, featol)
myprintf(                  '      d1max   = %8.1e', gamma)
myprintf(  '\nz0min    = %8g     opttol   = %8.1e', z0min, opttol)
myprintf(                  '      d2max   = %8.1e', delta)
myprintf(  '\nmu0      = %8.1e     steptol  = %8g', mu0  , steptol)
myprintf(                  '     bigcenter= %8g'  , bigcenter)

myprintf('\n\nLSQR:')
myprintf('\natol1    = %8.1e     atol2    = %8.1e', atol1 , atol2 )
myprintf(                  '      btol    = %8.1e', btol )
myprintf('\nconlim   = %8.1e     itnlim   = %8g'  , conlim, itnlim)
myprintf(                  '      show    = %8g'  , show )

% Method =  3;  %%% Hardwire LSQR
% Method = 41;  %%% Hardwire K4 and sparse LU

myprintf('\n\nMethod   = %8g     (1=chol  2=QR  3=LSQR  41=K4)', Method)

if wait
    myprintf('\n\nReview parameters... then type "return"\n')
    keyboard
end

if eta < 0
    myprintf('\n\nLinesearch disabled by eta < 0')
end

%---------------------------------------------------------------------
% All parameters have now been set.
% Check for valid Method.
%---------------------------------------------------------------------
time    = cputime;

if operator
    if Method==3
        % relax
    else
        myprintf('\n\nWhen A is an operator, we have to use Method = 3')
        Method = 3;
    end
end

if     Method== 1,  solver  = '  Chol';  head3 = '     Chol';
elseif Method== 2,  solver  = '    QR';  head3 = '       QR';
elseif Method== 3,  solver  = '  LSQR';  head3 = '  atol   LSQR Inexact';
elseif Method==41,  solver  = '    LU';  head3 = '        L       U   res';
else   error('Method must be 1, 2, 3, or 41')
end

%---------------------------------------------------------------------
% Categorize bounds and allow for fixed variables by modifying b.
%---------------------------------------------------------------------
[low,upp,fix] = pdxxxbounds( bl,bu );

nfix = length(fix);
if nfix > 0
    x1 = zn;   x1(fix) = bl(fix);
    r1 = pdxxxmat( Aname, 1, m, n, x1 );
    b  = b - r1;
    % At some stage, might want to look at normfix = norm(r1,inf);
end

%---------------------------------------------------------------------
% Scale the input data.
% The scaled variables are
%    xbar     = x/beta,
%    ybar     = y/zeta,
%    zbar     = z/zeta.
% Define
%    theta    = beta*zeta;
% The scaled function is
%    phibar   = ( 1   /theta) fbar(beta*xbar),
%    gradient = (beta /theta) grad,
%    Hessian  = (beta2/theta) hess.
%---------------------------------------------------------------------
beta   = xsize;   if beta==0, beta = 1; end    % beta scales b, x.
zeta   = zsize;   if zeta==0, zeta = 1; end    % zeta scales y, z.
theta  = beta*zeta;                            % theta scales obj.
% (theta could be anything, but theta = beta*zeta makes
% scaled grad = grad/zeta = 1 approximately if zeta is chosen right.)

bl(fix)= bl(fix)/beta;
bu(fix)= bu(fix)/beta;
bl(low)= bl(low)/beta;
bu(upp)= bu(upp)/beta;
d1     = d1*( beta/sqrt(theta) );
d2     = d2*( sqrt(theta)/beta );

beta2  = beta^2;
b      = b /beta;   y0     = y0/zeta;
x0     = x0/beta;   z0     = z0/zeta;

%---------------------------------------------------------------------
% Initialize vectors that are not fully used if bounds are missing.
%---------------------------------------------------------------------
rL     = zn;   rU    = zn;
cL     = zn;   cU    = zn;
x1     = zn;   x2    = zn;
z1     = zn;   z2    = zn;
dx1    = zn;   dx2   = zn;
dz1    = zn;   dz2   = zn;
clear zn

%---------------------------------------------------------------------
% Initialize x, y, z1, z2, objective, etc.
% 10 Aug 2003: z isn't needed here -- just at end for output.
%---------------------------------------------------------------------
x      = x0;
y      = y0;
x(fix) = bl(fix);
x(low) = max( x(low)          , bl(low));
x(upp) = min( x(upp)          , bu(upp));
x1(low)= max( x(low) - bl(low), x0min  );
x2(upp)= max(bu(upp) -  x(upp), x0min  );
z1(low)= max( z0(low)         , z0min  );
z2(upp)= max(-z0(upp)         , z0min  );
clear x0 y0 z0

%%%%%%%%%% Assume hess is diagonal for now. %%%%%%%%%%%%%%%%%%

if explicitF
    obj = (Fname'*x)*beta;  grad = Fname;   hess = zeros(n,1);
else
    [obj,grad,hess] = feval( Fname, (x*beta) );
end
obj  = obj        /theta;              % Scaled obj.
grad = grad*(beta /theta) + (d1.^2).*x;% grad includes x regularization.
H    = hess*(beta2/theta) + (d1.^2);   % H    includes x regularization.

%---------------------------------------------------------------------
% Compute primal and dual residuals:
%    r1 =  b - A*x - d2.^2*y
%    r2 =  grad - A'*y + (z2-z1)
%    rL =  bl - x + x1
%    rU = -bu + x + x2
%---------------------------------------------------------------------
[r1,r2,rL,rU,Pinf,Dinf] = ...
    pdxxxresid1( Aname,fix,low,upp, ...
    b,bl,bu,d1,d2,grad,rL,rU,x,x1,x2,y,z1,z2 );

%---------------------------------------------------------------------
% Initialize mu and complementarity residuals:
%    cL   = mu*e - X1*z1.
%    cU   = mu*e - X2*z2.
%
% 25 Jan 2001: Now that b and obj are scaled (and hence x,y,z),
%              we should be able to use mufirst = mu0 (absolute value).
%              0.1 worked poorly on StarTest1 with x0min = z0min = 0.1.
% 29 Jan 2001: We might as well use mu0 = x0min * z0min;
%              so that most variables are centered after a warm start.
% 29 Sep 2002: Use mufirst = mu0*(x0min * z0min),
%              regarding mu0 as a scaling of the initial center.
% 07 Aug 2003: mulast is controlled by opttol.
%              mufirst should not be smaller.
% 10 Aug 2003: Revert to mufirst = mu0 (absolute value).
%---------------------------------------------------------------------
% mufirst = mu0*(x0min * z0min);
mufirst = mu0;
mulast  = 0.1 * opttol;
mufirst = max( mufirst, mulast );
mu      = mufirst;
[cL,cU,center,Cinf,Cinf0] = ...
    pdxxxresid2( mu,low,upp,cL,cU,x1,x2,z1,z2 );
fmerit = pdxxxmerit( low,upp,r1,r2,rL,rU,cL,cU );

% Initialize other things.

PDitns    = 0;
converged = 0;
atol      = atol1;
atol2     = max( atol2, atolmin );
atolmin   = atol2;
pdDDD2    = d2;    % Global vector for diagonal matrix D2

%  Iteration log.

stepx   = 0;
stepz   = 0;
nf      = 0;
itncg   = 0;
nfail   = 0;

head1   = '\n\nItn   mu stepx stepz  Pinf  Dinf';
head2   = '  Cinf   Objective    nf  center';
myprintf([ head1 head2 head3 ])

regterm = norm(d1.*x)^2  +  norm(d2.*y)^2;
objreg  = obj  +  0.5*regterm;
objtrue = objreg * theta;

myprintf('\n%3g                 ', PDitns        )
myprintf('%6.1f%6.1f' , log10(Pinf ), log10(Dinf))
myprintf('%6.1f%15.7e', log10(Cinf0), objtrue    )
myprintf('   %8.1f'   , center                   )
if Method==41
    myprintf('   thresh=%7.1e', thresh)
end

if kminor
    myprintf('\n\nStart of first minor itn...\n')
    keyboard
end

%---------------------------------------------------------------------
% Main loop.
%---------------------------------------------------------------------

while ~converged
    PDitns = PDitns + 1;
    
    % 31 Jan 2001: Set atol according to progress, a la Inexact Newton.
    % 07 Feb 2001: 0.1 not small enough for Satellite problem.  Try 0.01.
    % 25 Apr 2001: 0.01 seems wasteful for Star problem.
    %              Now that starting conditions are better, go back to 0.1.
    
    r3norm = max([Pinf  Dinf  Cinf]);
    atol   = min([atol  r3norm*0.1]);
    atol   = max([atol  atolmin   ]);
    
    if Method<=3
        %-----------------------------------------------------------------
        %  Solve (*) for dy.
        %-----------------------------------------------------------------
        %  Define a damped Newton iteration for solving f = 0,
        %  keeping  x1, x2, z1, z2 > 0.  We eliminate dx1, dx2, dz1, dz2
        %  to obtain the system
        %
        %  [-H2  A' ] [dx] = [w ],   H2 = H + D1^2 + X1inv Z1 + X2inv Z2,
        %  [ A  D2^2] [dy] = [r1]    w  = r2 - X1inv(cL + Z1 rL)
        %                                    + X2inv(cU + Z2 rU),
        %
        %  which is equivalent to the least-squares problem
        %
        %     min || [ D A']dy  -  [  D w   ] ||,   D = H2^{-1/2}.     (*)
        %         || [  D2 ]       [D2inv r1] ||
        %-----------------------------------------------------------------
        H(low)  = H(low) + z1(low)./x1(low);
        H(upp)  = H(upp) + z2(upp)./x2(upp);
        w       = r2;
        w(low)  = w(low) - (cL(low) + z1(low).*rL(low))./x1(low);
        w(upp)  = w(upp) + (cU(upp) + z2(upp).*rU(upp))./x2(upp);
        
        H      = 1./H;    % H is now Hinv (NOTE!)
        H(fix) = 0;
        D      = sqrt(H);
        pdDDD1 = D;
        rw     = [explicitA Method m n 0 0 0];  % Passed to LSQR.
        
        if Method==1
            % --------------------------------------------------------------
            % Use chol to get dy
            % --------------------------------------------------------------
            AD   = Aname*sparse( 1:n, 1:n, D, n, n );
            ADDA = AD*AD' + sparse( 1:m, 1:m, (d2.^2), m, m );
            if PDitns==1, P = symamd(ADDA); end % Do ordering only once.
            
            [R,indef] = chol(ADDA(P,P));
            if indef
                myprintf('\n\n   chol says AD^2A'' is not pos def')
                myprintf('\n   Use bigger d2, or set options.Method = 2 or 3')
                break
            end
            
            % dy = ADDA \ rhs;
            rhs   = Aname*(H.*w) + r1;
            dy    = R \ (R'\rhs(P));
            dy(P) = dy;
            
            
        elseif Method==2
            % --------------------------------------------------------------
            % Use QR to get dy
            % --------------------------------------------------------------
            DAt  = sparse( 1:n, 1:n, D, n, n )*(Aname');
            if PDitns==1, P = colamd(DAt); end % Do ordering only once.
            
            if length(d2)==1
                D2 = d2*speye(m);
            else
                D2 = spdiags(d2,0,m,m);
            end
            
            DAt  = [ DAt;  D2     ];
            rhs  = [ D.*w; r1./d2 ];
            
            % dy = DAt \ rhs;
            [rhs,R] = qr(DAt(:,P),rhs,0);
            dy      = R \ rhs;
            dy(P)   = dy;
            
            
        else
            % --------------------------------------------------------------
            % Method=3.  Use LSQR (iterative solve) to get dy
            % --------------------------------------------------------------
            rhs    = [ D.*w; r1./d2 ];
            damp   = 0;
            
            if explicitA            % A is a sparse matrix.
                precon   = true;
                if precon    % Construct diagonal preconditioner for LSQR
                    AD     = Aname*sparse( 1:n, 1:n, D, n, n );
                    AD     = AD.^2;
                    wD     = sum(AD,2); % Sum of sqrs of each row of AD.  %(Sparse)
                    wD     = sqrt( full(wD) + (d2.^2) );                  %(Dense)
                    pdDDD3 = 1./wD;
                    clear AD wD
                end
            else                    % A is an operator
                precon   = false;
            end
            
            rw(7)        = precon;
            info.atolmin = atolmin;
            info.r3norm  = fmerit;  % Must be the 2-norm here.
            
            [ dy, istop, itncg, outfo ] = ...
                pdxxxlsqr( nb,m,'pdxxxlsqrmat',Aname,rw,rhs,damp, ...
                atol,btol,conlim,itnlim,show,info );
            
            if precon, dy = pdDDD3 .* dy; end
            
            if istop==3  |  istop==7   % conlim or itnlim
                myprintf('\n    LSQR stopped early:  istop = %3d', istop)
            end
            
            atolold   = outfo.atolold;
            atol      = outfo.atolnew;
            r3ratio   = outfo.r3ratio;
            CGitns    = CGitns + itncg;
        end % computation of dy
        
        % dy is now known.  Get dx, dx1, dx2, dz1, dz2.
        
        grad      =   pdxxxmat( Aname,2,m,n,dy );  % grad = A'dy
        grad(fix) =   0;    % Is this needed?      % grad is a work vector
        dx        =   H .* (grad - w);
        dx1(low)  = - rL(low) + dx(low);
        dx2(upp)  = - rU(upp) - dx(upp);
        dz1(low)  =  (cL(low) - z1(low).*dx1(low)) ./ x1(low);
        dz2(upp)  =  (cU(upp) - z2(upp).*dx2(upp)) ./ x2(upp);
        
        
    elseif Method==41
        %-----------------------------------------------------------------
        % Solve the symmetric-structure 4 x 4 system K4 dv = t:
        %
        %     ( X1     Z1      ) [dz1] = [tL],    tL = cL + Z1 rL
        %     (    X2 -Z2      ) [dz2]   [tU]     tU = cU + Z2 rU
        %     ( I  -I -H1  A'  ) [dx ]   [r2]
        %     (        A  D2^2 ) [dy ]   [r1]
        %-----------------------------------------------------------------
        X1  = ones(n,1);   X1(low) = x1(low);
        X2  = ones(n,1);   X2(upp) = x2(upp);
        if length(d2)==1
            D22 = d2^2*speye(m);
        else
            D22 = spdiags(d2.^2,0,m,m);
        end
        
        Onn = sparse(n,n);
        Omn = sparse(m,n);
        
        if PDitns==1    % First time through: Choose LU ordering from
            % lower-triangular part of dummy K4
            I1 = sparse( low, low, ones(length(low),1), n, n );
            I2 = sparse( upp, upp, ones(length(upp),1), n, n );
            K4 =[speye(n)            Onn             Onn                Omn'
                Onn          speye(n)           Onn                Omn'
                I1              I2              speye(n)           Omn'
                Omn             Omn             Aname         speye(m)];
            p      = symamd(K4);
            % mydisp(' ');  keyboard
        end
        
        K4  = [spdiags(X1,0,n,n)   Onn             spdiags(z1,0,n,n)  Omn'
            Onn          spdiags(X2,0,n,n) -spdiags(z2,0,n,n)  Omn'
            I1             -I2             -spdiags(H,0,n,n)   Aname'
            Omn             Omn             Aname              D22];
        
        tL  = zeros(n,1);   tL(low) = cL(low) + z1(low).*rL(low);
        tU  = zeros(n,1);   tU(upp) = cU(upp) + z2(upp).*rU(upp);
        rhs = [tL; tU; r2; r1];
        
        %   dv  = K4 \ rhs;        % BIG SYSTEM!
        
        [L,U,P] = lu( K4(p,p), thresh );     % P K4 = L U
        dv      = U \ (L \ (P*rhs(p)));
        dv(p)   = dv;
        resK4   = norm((K4*dv - rhs),inf) / norm(rhs,inf);
        
        dz1 = dv(    1:  n);
        dz2 = dv(  n+1:2*n);
        dx  = dv(2*n+1:3*n);
        dy  = dv(3*n+1:3*n+m);
        
        dx1(low) = - rL(low) + dx(low);
        dx2(upp) = - rU(upp) - dx(upp);
    end
    
    
    %-------------------------------------------------------------------
    % Find the maximum step.
    % 13 Aug 2003: We need stepxL, stepxU also to keep x feasible
    %              so that nonlinear functions are defined.
    % 18 Nov 2003: But this gives stepx = 0 for lptest.  (??)
    %--------------------------------------------------------------------
    stepx1 = pdxxxstep( x1(low), dx1(low) );
    stepx2 = pdxxxstep( x2(upp), dx2(upp) );
    stepz1 = pdxxxstep( z1(low), dz1(low) );
    stepz2 = pdxxxstep( z2(upp), dz2(upp) );
    %  stepxL = pdxxxstep(  x(low),  dx(low) );
    %  stepxU = pdxxxstep(  x(upp),  dx(upp) );
    %  stepx  = min( [stepx1,   stepx2,   stepxL,   stepxU] );
    stepx  = min( [stepx1,   stepx2] );
    stepz  = min( [stepz1,   stepz2] );
    stepx  = min( [steptol*stepx, 1] );
    stepz  = min( [steptol*stepz, 1] );
    if stepSame                      % For NLPs, force same step
        stepx = min( stepx, stepz );   % (true Newton method)
        stepz = stepx;
    end
    
    %-------------------------------------------------------------------
    % Backtracking linesearch.
    %-------------------------------------------------------------------
    fail     =  true;
    nf       =  0;
    
    while nf < maxf
        nf      = nf + 1;
        x       = x        +  stepx * dx;
        y       = y        +  stepz * dy;
        x1(low) = x1(low)  +  stepx * dx1(low);
        x2(upp) = x2(upp)  +  stepx * dx2(upp);
        z1(low) = z1(low)  +  stepz * dz1(low);
        z2(upp) = z2(upp)  +  stepz * dz2(upp);
        if explicitF
            obj = (Fname'*x)*beta;  grad = Fname;   hess = zeros(n,1);
        else
            [obj,grad,hess] = feval( Fname, (x*beta) );
        end
        obj        = obj /theta;
        grad       = grad*(beta /theta)  +  (d1.^2).*x;
        H          = hess*(beta2/theta)  +  (d1.^2);
        
        [r1,r2,rL,rU,Pinf,Dinf] = ...
            pdxxxresid1( Aname,fix,low,upp, ...
            b,bl,bu,d1,d2,grad,rL,rU,x,x1,x2,y,z1,z2 );
        
        [cL,cU,center,Cinf,Cinf0] = ...
            pdxxxresid2( mu,low,upp,cL,cU,x1,x2,z1,z2 );
        fmeritnew = pdxxxmerit( low,upp,r1,r2,rL,rU,cL,cU );
        
        step      = min( stepx, stepz );
        
        if fmeritnew <= (1 - eta*step)*fmerit
            fail = false;
            break;
        end
        
        % Merit function didn't decrease.
        % Restore variables to previous values.
        % (This introduces a little error, but save lots of space.)
        
        x       = x        -  stepx * dx;
        y       = y        -  stepz * dy;
        x1(low) = x1(low)  -  stepx * dx1(low);
        x2(upp) = x2(upp)  -  stepx * dx2(upp);
        z1(low) = z1(low)  -  stepz * dz1(low);
        z2(upp) = z2(upp)  -  stepz * dz2(upp);
        
        % Back-track.
        % If it's the first time,
        % make stepx and stepz the same.
        
        if nf==1 & stepx~=stepz
            stepx = step;
        elseif nf < maxf
            stepx = stepx/2;
        end;
        stepz = stepx;
    end
    
    if fail
        myprintf('\n     Linesearch failed (nf too big)');
        nfail = nfail + 1;
    else
        nfail = 0;
    end
    
    %-------------------------------------------------------------------
    % Set convergence measures.
    %--------------------------------------------------------------------
    regterm = norm(d1.*x)^2  +  norm(d2.*y)^2;
    objreg  = obj  +  0.5*regterm;
    objtrue = objreg * theta;
    
    primalfeas    = Pinf  <=  featol;
    dualfeas      = Dinf  <=  featol;
    complementary = Cinf0 <=  opttol;
    enough        = PDitns>=       4;  % Prevent premature termination.
    converged     = primalfeas  &  dualfeas  &  complementary  &  enough;
    
    %-------------------------------------------------------------------
    % Iteration log.
    %-------------------------------------------------------------------
    str1    = sprintf('\n%3g%5.1f' , PDitns      , log10(mu)   );
    str2    = sprintf('%6.3f%6.3f' , stepx       , stepz       );
    if stepx < 0.001 | stepz < 0.001
        str2 = sprintf('%6.1e%6.1e' , stepx       , stepz       );
    end
    
    str3    = sprintf('%6.1f%6.1f' , log10(Pinf) , log10(Dinf));
    str4    = sprintf('%6.1f%15.7e', log10(Cinf0), objtrue     );
    str5    = sprintf('%3g%8.1f'   , nf          , center      );
    if center > 99999
        str5 = sprintf('%3g%8.1e'   , nf          , center      );
    end
    myprintf([str1 str2 str3 str4 str5])
    
    if     Method== 1
        if PDitns==1, myprintf(' %8g', nnz(R)); end
    elseif Method== 2
        if PDitns==1, myprintf(' %8g', nnz(R)); end
    elseif Method== 3
        myprintf(' %5.1f%7g%7.3f', log10(atolold), itncg, r3ratio)
    elseif Method==41,
        resK4 = max( resK4, 1e-99 );
        myprintf(' %8g%8g%6.1f', nnz(L),nnz(U),log10(resK4))
    end
    
    %-------------------------------------------------------------------
    % Test for termination.
    %-------------------------------------------------------------------
    if kminor
        myprintf( '\nStart of next minor itn...\n')
        keyboard
    end
    
    if converged
        myprintf('\nConverged')
    elseif PDitns >= maxitn
        myprintf('\nToo many iterations')
        inform = 1;
        break
    elseif nfail  >= maxfail
        myprintf('\nToo many linesearch failures')
        inform = 2;
        break
    elseif step   <= 1e-10
        myprintf('\nStep lengths too small')
        inform = 3;
        break
    else
        
        % Reduce mu, and reset certain residuals.
        
        stepmu  = min( stepx , stepz   );
        stepmu  = min( stepmu, steptol );
        muold   = mu;
        mu      = mu   -  stepmu * mu;
        if center >= bigcenter,  mu = muold;  end
        
        % mutrad = mu0*(sum(Xz)/n); % 24 May 1998: Traditional value, but
        % mu     = min(mu,mutrad ); % it seemed to decrease mu too much.
        
        mu      = max(mu,mulast);  % 13 Jun 1998: No need for smaller mu.
        [cL,cU,center,Cinf,Cinf0] = ...
            pdxxxresid2( mu,low,upp,cL,cU,x1,x2,z1,z2 );
        fmerit = pdxxxmerit( low,upp,r1,r2,rL,rU,cL,cU );
        
        % Reduce atol for LSQR (and SYMMLQ).
        % NOW DONE AT TOP OF LOOP.
        
        atolold = atol;
        
        % if atol > atol2
        %   atolfac = (mu/mufirst)^0.25;
        %   atol    = max( atol*atolfac, atol2 );
        % end
        
        % atol = min( atol, mu );     % 22 Jan 2001: a la Inexact Newton.
        % atol = min( atol, 0.5*mu ); % 30 Jan 2001: A bit tighter
        
        % If the linesearch took more than one function (nf > 1),
        % we assume the search direction needed more accuracy
        % (though this may be true only for LPs).
        % 12 Jun 1998: Ask for more accuracy if nf > 2.
        % 24 Nov 2000: Also if the steps are small.
        % 30 Jan 2001: Small steps might be ok with warm start.
        % 06 Feb 2001: Not necessarily.  Reinstated tests in next line.
        
        if nf > 2  |  step <= 0.01
            atol = atolold*0.1;
        end
    end
end
%---------------------------------------------------------------------
% End of main loop.
%---------------------------------------------------------------------

% Print statistics.

x(fix) = 0;                 % Exclude x(fix) temporarily from |x|.
z      = zeros(n,1);        % Exclude z(fix) also.
z(low) = z1(low);
z(upp) = z(upp) - z2(upp);

myprintf('\n\nmax |x| =%10.3f', norm(x,inf))
myprintf('    max |y| =%10.3f', norm(y,inf))
myprintf('    max |z| =%10.3f', norm(z,inf))  % excludes z(fix)
myprintf('   scaled')

bl(fix)= bl(fix)*beta;      % Unscale bl, bu, x, y, z.
bu(fix)= bu(fix)*beta;
bl(low)= bl(low)*beta;
bu(upp)= bu(upp)*beta;

x = x*beta;   y = y*zeta;   z = z*zeta;

myprintf(  '\nmax |x| =%10.3f', norm(x,inf))
myprintf('    max |y| =%10.3f', norm(y,inf))
myprintf('    max |z| =%10.3f', norm(z,inf))  % excludes z(fix)
myprintf(' unscaled')

x(fix) = bl(fix);           % Reinstate x(fix).

% Reconstruct b.
b      = b *beta;
if nfix > 0
    x1 = zeros(n,1);   x1(fix) = bl(fix);
    r1 = pdxxxmat( Aname, 1, m, n, x1 );
    b  = b + r1;
    myprintf('\nmax |x| and max |z| exclude fixed variables')
end

% Evaluate function at final point.
% Reconstruct z.  This finally defines z(fix).
if explicitF
    obj = (Fname'*x);  grad = Fname;   hess = zeros(n,1);
else
    [obj,grad,hess] = feval( Fname, x );
end
z      = grad - pdxxxmat( Aname,2,m,n,y );  % z = grad - A'y

time   = cputime - time;
str1   = sprintf('\nPDitns  =%10g', PDitns );
str2   = sprintf(     'itns =%10g', CGitns );
myprintf( [str1 ' ' solver str2] )
myprintf('    time    =%10.1f', time);

pdxxxdistrib( abs(x),abs(z) );   % Private function

if wait
    keyboard
end
%-----------------------------------------------------------------------
% End function pdco.m
%-----------------------------------------------------------------------


function [low,upp,fix] = pdxxxbounds( bl,bu )

% Categorize various types of bounds.
% pos overlaps with low.
% neg overlaps with upp.
% two overlaps with low and upp.
% fix and free are disjoint from all other sets.

bigL = -9.9e+19;
bigU =  9.9e+19;
pos  =  find( bl==0    & bu>=bigU );
neg  =  find( bl<=bigL & bu==0    );
low  =  find( bl> bigL & bl< bu   );
upp  =  find( bu< bigU & bl< bu   );
two  =  find( bl> bigL & bu< bigU & bl< bu );
fix  =  find( bl==bu );
free =  find( bl<=bigL & bu>=bigU );

myprintf('\n\nBounds:\n  [0,inf]  [-inf,0]')
myprintf('  Finite bl  Finite bu  Two bnds   Fixed    Free')
myprintf('\n %8g %9g %10g %10g %9g %7g %7g',     ...
    length(pos), length(neg), length(low),  ...
    length(upp), length(two), length(fix), length(free))

%-----------------------------------------------------------------------
% End private function pdxxxbounds
%-----------------------------------------------------------------------


function pdxxxdistrib( x,z )

% pdxxxdistrib(x) or pdxxxdistrib(x,z) prints the
% distribution of 1 or 2 vectors.
%
% 18 Dec 2000.  First version with 2 vectors.

two  = nargin > 1;
myprintf('\n\nDistribution of vector     x')
if two, myprintf('         z'); end

x1   = 10^(floor(log10(max(x)+eps)) + 1);
z1   = 10^(floor(log10(max(z)+eps)) + 1);
x1   = max(x1,z1);
kmax = 10;

for k = 1:kmax
    x2 = x1;    x1 = x1/10;
    if k==kmax, x1 = 0; end
    nx = length(find(x>=x1 & x<x2));
    myprintf('\n[%7.3g,%7.3g )%10g', x1, x2, nx);
    if two
        nz = length(find(z>=x1 & z<x2));
        myprintf('%10g', nz);
    end
end

mydisp(' ')

%-----------------------------------------------------------------------
% End private function pdxxxdistrib
%-----------------------------------------------------------------------


function [ x, istop, itn, outfo ] = ...
    pdxxxlsqr( m, n, aprodname, iw, rw, b, damp, ...
    atol, btol, conlim, itnlim, show, info )

% Special version of LSQR for use with pdco.m.
% It continues with a reduced atol if a pdco-specific test isn't
% satisfied with the input atol.
%
% LSQR solves  Ax = b  or  min ||b - Ax||_2  if damp = 0,
% or   min || (b)  -  (  A   )x ||   otherwise.
%          || (0)     (damp I)  ||2
% A  is an m by n matrix defined by  y = aprod( mode,m,n,x,iw,rw ),
% where the parameter 'aprodname' refers to a function 'aprod' that
% performs the matrix-vector operations.
% If mode = 1,   aprod  must return  y = Ax   without altering x.
% If mode = 2,   aprod  must return  y = A'x  without altering x.
% WARNING:   The file containing the function 'aprod'
%            must not be called aprodname.m !!!!

%-----------------------------------------------------------------------
% LSQR uses an iterative (conjugate-gradient-like) method.
% For further information, see
% 1. C. C. Paige and M. A. Saunders (1982a).
%    LSQR: An algorithm for sparse linear equations and sparse least squares,
%    ACM TOMS 8(1), 43-71.
% 2. C. C. Paige and M. A. Saunders (1982b).
%    Algorithm 583.  LSQR: Sparse linear equations and least squares problems,
%    ACM TOMS 8(2), 195-209.
% 3. M. A. Saunders (1995).  Solution of sparse rectangular systems using
%    LSQR and CRAIG, BIT 35, 588-604.
%
% Input parameters:
% iw, rw      are not used by lsqr, but are passed to aprod.
% atol, btol  are stopping tolerances.  If both are 1.0e-9 (say),
%             the final residual norm should be accurate to about 9 digits.
%             (The final x will usually have fewer correct digits,
%             depending on cond(A) and the size of damp.)
% conlim      is also a stopping tolerance.  lsqr terminates if an estimate
%             of cond(A) exceeds conlim.  For compatible systems Ax = b,
%             conlim could be as large as 1.0e+12 (say).  For least-squares
%             problems, conlim should be less than 1.0e+8.
%             Maximum precision can be obtained by setting
%             atol = btol = conlim = zero, but the number of iterations
%             may then be excessive.
% itnlim      is an explicit limit on iterations (for safety).
% show = 1    gives an iteration log,
% show = 0    suppresses output.
% info        is a structure special to pdco.m, used to test if
%             was small enough, and continuing if necessary with smaller atol.
%
%
% Output parameters:
% x           is the final solution.
% istop       gives the reason for termination.
% istop       = 1 means x is an approximate solution to Ax = b.
%             = 2 means x approximately solves the least-squares problem.
% rnorm       = norm(r) if damp = 0, where r = b - Ax,
%             = sqrt( norm(r)**2  +  damp**2 * norm(x)**2 ) otherwise.
% xnorm       = norm(x).
% var         estimates diag( inv(A'A) ).  Omitted in this special version.
% outfo       is a structure special to pdco.m, returning information
%             about whether atol had to be reduced.
%
% Other potential output parameters:
% anorm, acond, arnorm, xnorm
%
%        1990: Derived from Fortran 77 version of LSQR.
% 22 May 1992: bbnorm was used incorrectly.  Replaced by anorm.
% 26 Oct 1992: More input and output parameters added.
% 01 Sep 1994: Matrix-vector routine is now a parameter 'aprodname'.
%              Print log reformatted.
% 14 Jun 1997: show  added to allow printing or not.
% 30 Jun 1997: var   added as an optional output parameter.
%              It returns an estimate of diag( inv(A'A) ).
% 12 Feb 2001: atol  can now be reduced and iterations continued if necessary.
%              info, outfo are new problem-dependent parameters for such purposes.
%              In this version they are specialized for pdco.m.
%-----------------------------------------------------------------------

%     Initialize.

msg=['The exact solution is  x = 0                              '
    'Ax - b is small enough, given atol, btol                  '
    'The least-squares solution is good enough, given atol     '
    'The estimate of cond(Abar) has exceeded conlim            '
    'Ax - b is small enough for this machine                   '
    'The least-squares solution is good enough for this machine'
    'Cond(Abar) seems to be too large for this machine         '
    'The iteration limit has been reached                      '];

% wantvar= nargout >= 6;
% if wantvar, var = zeros(n,1); end

itn    = 0;           istop  = 0;             nstop  = 0;
ctol   = 0;           if conlim > 0, ctol = 1/conlim; end;
anorm  = 0;           acond  = 0;
dampsq = damp^2;      ddnorm = 0;             res2   = 0;
xnorm  = 0;           xxnorm = 0;             z      = 0;
cs2    = -1;          sn2    = 0;

% Set up the first vectors u and v for the bidiagonalization.
% These satisfy  beta*u = b,  alfa*v = A'u.

u      = b(1:m);      x    = zeros(n,1);
alfa   = 0;           beta = norm( u );
if beta > 0
    u = (1/beta) * u;   v = feval( aprodname, 2, m, n, u, iw, rw );
    alfa = norm( v );
end
if alfa > 0
    v = (1/alfa) * v;    w = v;
end

arnorm = alfa * beta;   if arnorm==0, mydisp(msg(1,:)); return, end

rhobar = alfa;                phibar = beta;          bnorm  = beta;
rnorm  = beta;
head1  = '   Itn      x(1)      Function';
head2  = ' Compatible   LS      Norm A   Cond A';

if show
    mydisp(' ')
    mydisp([head1 head2])
    test1  = 1;         test2  = alfa / beta;
    str1   = sprintf( '%6g %12.5e %10.3e',   itn, x(1), rnorm );
    str2   = sprintf( '  %8.1e %8.1e',       test1, test2 );
    mydisp([str1 str2])
end

%----------------------------------------------------------------
% Main iteration loop.
%----------------------------------------------------------------
while itn < itnlim
    itn = itn + 1;
    % Perform the next step of the bidiagonalization to obtain the
    % next beta, u, alfa, v.  These satisfy the relations
    % beta*u  =  A*v   -  alfa*u,
    % alfa*v  =  A'*u  -  beta*v.
    
    u    = feval( aprodname, 1, m, n, v, iw, rw )  -  alfa*u;
    beta = norm( u );
    if beta > 0
        u     = (1/beta) * u;
        anorm = norm([anorm alfa beta damp]);
        v     = feval( aprodname, 2, m, n, u, iw, rw )  -  beta*v;
        alfa  = norm( v );
        if alfa > 0,  v = (1/alfa) * v; end
    end
    
    % Use a plane rotation to eliminate the damping parameter.
    % This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
    
    rhobar1 = norm([rhobar damp]);
    cs1     = rhobar / rhobar1;
    sn1     = damp   / rhobar1;
    psi     = sn1 * phibar;
    phibar  = cs1 * phibar;
    
    % Use a plane rotation to eliminate the subdiagonal element (beta)
    % of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
    
    rho     = norm([rhobar1 beta]);
    cs      =   rhobar1/ rho;
    sn      =   beta   / rho;
    theta   =   sn * alfa;
    rhobar  = - cs * alfa;
    phi     =   cs * phibar;
    phibar  =   sn * phibar;
    tau     =   sn * phi;
    
    % Update x and w.
    
    t1      =   phi  /rho;
    t2      = - theta/rho;
    dk      =   (1/rho)*w;
    
    x       = x      +  t1*w;
    w       = v      +  t2*w;
    ddnorm  = ddnorm +  norm(dk)^2;
    % if wantvar, var = var  +  dk.*dk; end
    
    % Use a plane rotation on the right to eliminate the
    % super-diagonal element (theta) of the upper-bidiagonal matrix.
    % Then use the result to estimate  norm(x).
    
    delta   =   sn2 * rho;
    gambar  = - cs2 * rho;
    rhs     =   phi  -  delta * z;
    zbar    =   rhs / gambar;
    xnorm   =   sqrt(xxnorm + zbar^2);
    gamma   =   norm([gambar theta]);
    cs2     =   gambar / gamma;
    sn2     =   theta  / gamma;
    z       =   rhs    / gamma;
    xxnorm  =   xxnorm  +  z^2;
    
    % Test for convergence.
    % First, estimate the condition of the matrix  Abar,
    % and the norms of  rbar  and  Abar'rbar.
    
    acond   =   anorm * sqrt( ddnorm );
    res1    =   phibar^2;
    res2    =   res2  +  psi^2;
    rnorm   =   sqrt( res1 + res2 );
    arnorm  =   alfa * abs( tau );
    
    % Now use these norms to estimate certain other quantities,
    % some of which will be small near a solution.
    
    test1   =   rnorm / bnorm;
    test2   =   arnorm/( anorm * rnorm );
    test3   =       1 / acond;
    t1      =   test1 / (1    +  anorm * xnorm / bnorm);
    rtol    =   btol  +  atol *  anorm * xnorm / bnorm;
    
    % The following tests guard against extremely small values of
    % atol, btol  or  ctol.  (The user may have set any or all of
    % the parameters  atol, btol, conlim  to 0.)
    % The effect is equivalent to the normal tests using
    % atol = eps,  btol = eps,  conlim = 1/eps.
    
    if itn >= itnlim,   istop = 7; end
    if 1 + test3  <= 1, istop = 6; end
    if 1 + test2  <= 1, istop = 5; end
    if 1 + t1     <= 1, istop = 4; end
    
    % Allow for tolerances set by the user.
    
    if  test3 <= ctol,  istop = 3; end
    if  test2 <= atol,  istop = 2; end
    if  test1 <= rtol,  istop = 1; end
    
    %-------------------------------------------------------------------
    % SPECIAL TEST THAT DEPENDS ON pdco.m.
    % Aname in pdco   is  iw in lsqr.
    % dy              is  x
    % Other stuff     is in info.
    % We allow for diagonal preconditioning in pdDDD3.
    %-------------------------------------------------------------------
    if istop > 0
        r3new     = arnorm;
        r3ratio   = r3new / info.r3norm;
        atolold   = atol;
        atolnew   = atol;
        
        if atol > info.atolmin
            if     r3ratio <= 0.1     % dy seems good
                % Relax
            elseif r3ratio <= 0.5     % Accept dy but make next one more accurate.
                atolnew = atolnew * 0.1;
            else                      % Recompute dy more accurately
                myprintf('\n                                ')
                myprintf('                                ')
                myprintf(' %5.1f%7g%7.3f', log10(atolold), itn, r3ratio)
                atol    = atol * 0.1;
                atolnew = atol;
                istop   = 0;
            end
        end
        
        outfo.atolold = atolold;
        outfo.atolnew = atolnew;
        outfo.r3ratio = r3ratio;
    end
    
    %-------------------------------------------------------------------
    % See if it is time to print something.
    %-------------------------------------------------------------------
    prnt = 0;
    if n     <= 40       , prnt = 1; end
    if itn   <= 10       , prnt = 1; end
    if itn   >= itnlim-10, prnt = 1; end
    if rem(itn,10)==0    , prnt = 1; end
    if test3 <=  2*ctol  , prnt = 1; end
    if test2 <= 10*atol  , prnt = 1; end
    if test1 <= 10*rtol  , prnt = 1; end
    if istop ~=  0       , prnt = 1; end
    
    if prnt==1
        if show
            str1 = sprintf( '%6g %12.5e %10.3e',   itn, x(1), rnorm );
            str2 = sprintf( '  %8.1e %8.1e',       test1, test2 );
            str3 = sprintf( ' %8.1e %8.1e',        anorm, acond );
            mydisp([str1 str2 str3])
        end
    end
    if istop > 0, break, end
end

% End of iteration loop.
% Print the stopping condition.

if show
    mydisp(' ')
    mydisp('LSQR finished')
    mydisp(msg(istop+1,:))
    mydisp(' ')
    str1 = sprintf( 'istop  =%8g   itn    =%8g',      istop, itn    );
    str2 = sprintf( 'anorm  =%8.1e   acond  =%8.1e',  anorm, acond  );
    str3 = sprintf( 'rnorm  =%8.1e   arnorm =%8.1e',  rnorm, arnorm );
    str4 = sprintf( 'bnorm  =%8.1e   xnorm  =%8.1e',  bnorm, xnorm  );
    mydisp([str1 '   ' str2])
    mydisp([str3 '   ' str4])
    mydisp(' ')
end

%-----------------------------------------------------------------------
% End private function pdxxxlsqr
%-----------------------------------------------------------------------


function y = pdxxxlsqrmat( mode, mlsqr, nlsqr, x, Aname, rw )

% pdxxxlsqrmat is required by pdco.m (when it calls pdxxxlsqr.m).
% It forms Mx or M'x for some operator M that depends on Method below.
%
% mlsqr, nlsqr  are the dimensions of the LS problem that lsqr is solving.
%
% Aname is pdco's Aname.
%
% rw contains parameters [explicitA Method LSdamp]
% from pdco.m to say which least-squares subproblem is being solved.
%
% global pdDDD1 pdDDD3 provides various diagonal matrices
% for each value of Method.

%-----------------------------------------------------------------------
% 17 Mar 1998: First version to go with pdsco.m and lsqr.m.
% 01 Apr 1998: global pdDDD1 pdDDD3 now used for efficiency.
% 11 Feb 2000: Added diagonal preconditioning for LSQR, solving for dy.
% 14 Dec 2000: Added diagonal preconditioning for LSQR, solving for dx.
% 12 Feb 2001: Included in pdco.m as private function.
%              Specialized to allow only solving for dy.
% 03 Oct 2002: First version to go with pdco.m with general H2 and D2.
% 16 Oct 2002: Aname is now the user's Aname.
%-----------------------------------------------------------------------

global pdDDD1 pdDDD2 pdDDD3

Method    = rw(2);
precon    = rw(7);

if Method==3
    % The operator M is [ D1 A'; D2 ].
    m = nlsqr;
    n = mlsqr - m;
    if mode==1
        if precon, x = pdDDD3.*x; end
        t = pdxxxmat( Aname, 2, m, n, x );   % Ask 'aprod' to form t = A'x.
        y = [ (pdDDD1.*t); (pdDDD2.*x) ];
    else
        t = pdDDD1.*x(1:n);
        y = pdxxxmat( Aname, 1, m, n, t );   % Ask 'aprod' to form y = A t.
        y = y + pdDDD2.*x(n+1:mlsqr);
        if precon, y = pdDDD3.*y; end
    end
    
else
    error('Error in pdxxxlsqrmat: Only Method = 3 is allowed at present')
end

%-----------------------------------------------------------------------
% End private function pdxxxlsqrmat
%-----------------------------------------------------------------------


function y = pdxxxmat( Aname, mode, m, n, x )

%        y = pdxxxmat( Aname, mode, m, n, x )
%    computes y = Ax (mode=1) or A'x (mode=2)
%    for a matrix A defined by pdco's input parameter Aname.

%-----------------------------------------------------------------------
% 04 Apr 1998: Default A*x and A'*y function for pdco.m.
%              Assumed A was a global matrix pdAAA created by pdco.m
%              from the user's input parameter A.
% 16 Oct 2002: pdAAA eliminated to save storage.
%              User's parameter Aname is now passed thru to here.
% 01 Nov 2002: Bug: feval had one too many parameters.
%-----------------------------------------------------------------------

if (ischar(Aname) || isa(Aname, 'function_handle'))
    y = feval( Aname, mode, m, n, x );
else
    if mode==1,  y = Aname*x;  else  y = Aname'*x;  end
end

%-----------------------------------------------------------------------
% End private function pdxxxmat
%-----------------------------------------------------------------------


function fmerit = pdxxxmerit( low,upp,r1,r2,rL,rU,cL,cU )

% Evaluate the merit function for Newton's method.
% It is the 2-norm of the three sets of residuals.

f = [norm(r1)
    norm(r2)
    norm(rL(low))
    norm(rU(upp))
    norm(cL(low))
    norm(cU(upp))];
fmerit = norm(f);

%-----------------------------------------------------------------------
% End private function pdxxxmerit
%-----------------------------------------------------------------------


function [r1,r2,rL,rU,Pinf,Dinf] =    ...
    pdxxxresid1( Aname,fix,low,upp, ...
    b,bl,bu,d1,d2,grad,rL,rU,x,x1,x2,y,z1,z2 )

% Form residuals for the primal and dual equations.
% rL, rU are output, but we input them as full vectors
% initialized (permanently) with any relevant zeros.
% 13 Aug 2003: z2-z1 coded more carefully
%              (although MATLAB was doing the right thing).
% 19 Nov 2003: r2(fix) = 0 has to be done after r2 = grad - r2;

m       = length(b);
n       = length(bl);
x(fix)  = 0;
r1      = pdxxxmat( Aname, 1, m, n, x );
r2      = pdxxxmat( Aname, 2, m, n, y );

r1      = b    - r1 - (d2.^2).*y;
r2      = grad - r2;  % + (z2-z1);        % grad includes (d1.^2)*x
r2(fix) = 0;
r2(upp) = r2(upp) + z2(upp);
r2(low) = r2(low) - z1(low);
rL(low) = (  bl(low) - x(low)) + x1(low);
rU(upp) = (- bu(upp) + x(upp)) + x2(upp);

Pinf    = max([norm(r1,inf) norm(rL(low),inf) norm(rU(upp),inf)]);
Dinf    =      norm(r2,inf);
Pinf    = max( Pinf, 1e-99 );
Dinf    = max( Dinf, 1e-99 );

%-----------------------------------------------------------------------
% End private function pdxxxresid1
%-----------------------------------------------------------------------


function [cL,cU,center,Cinf,Cinf0] = ...
    pdxxxresid2( mu,low,upp,cL,cU,x1,x2,z1,z2 )

% Form residuals for the complementarity equations.
% cL, cU are output, but we input them as full vectors
% initialized (permanently) with any relevant zeros.
% Cinf  is the complementarity residual for X1 z1 = mu e, etc.
% Cinf0 is the same for mu=0 (i.e., for the original problem).

x1z1    = x1(low).*z1(low);
x2z2    = x2(upp).*z2(upp);
cL(low) = mu - x1z1;
cU(upp) = mu - x2z2;

maxXz   = max( [max(x1z1) max(x2z2)] );
minXz   = min( [min(x1z1) min(x2z2)] );
maxXz   = max( maxXz, 1e-99 );
minXz   = max( minXz, 1e-99 );
center  = maxXz / minXz;
Cinf    = max([norm(cL(low),inf) norm(cU(upp),inf)]);
Cinf0   = maxXz;

%-----------------------------------------------------------------------
% End private function pdxxxresid2
%-----------------------------------------------------------------------


function step = pdxxxstep( x, dx )

% Assumes x > 0.
% Finds the maximum step such that x + step*dx >= 0.

step     = 1e+20;
blocking = find( dx < 0 );
if length( blocking ) > 0
    steps  = x(blocking) ./ (- dx(blocking));
    step   = min( steps );
end

%
% Copyright (c) 2006. Michael Saunders
%

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%

function myprintf(s,a,b,c,d,e,f,g,h,i,j,k,l,m,n)


function mydisp(s)

function options = pdcoSet(varargin)
%pdcoSet creates or alters an options structure for pdco.m.
%It is modeled after MATLAB's original version of optimset.
%
%   options = pdcoSet (with no input arguments) creates a structure
%   with all fields set to their default values.  Each field is an
%   option (also called a parameter).
%
%   pdcoSet (with no input or output arguments) displays all options
%   and their default values.
%
%   options = pdcoSet('PARAM1',VALUE1,'PARAM2',VALUE2,...);
%   creates a structure with the specified named parameters and values.
%   Any unspecified parameters are set to [].
%   It is sufficient to use only the leading characters that uniquely
%   identify the parameter.  Case is ignored for parameter names.
%
%   NOTE: For values that are strings, correct case and the complete string
%   are required; if an invalid string is provided, the default is used.
%
%   options = pdcoSet(oldopts,'PARAM1',VALUE1,...);
%   creates a copy of oldopts with the named parameters reset to the
%   specified values.
%
%   options = pdcoSet(oldopts,newopts);
%   combines an existing structure oldopts with a new structure newopts.
%   Any parameters in newopts with non-empty values overwrite the
%   corresponding old parameters in oldopts.
%
% options.MaxIter     Maximum iterations of the primal-dual barrier method.
%                     Most problems should solve within 30 PDitns.
% options.FeaTol      Accuracy for satisfying Ax + D2 r = b, A'y + z = gobj
%                     and x - x1 = bl, x + x2 = bu, where x1, x2 > 0.
%                     1e-6 is typically small enough.
%                     1e-5 may be acceptable also.
% options.OptTol      Accuracy for satisfying x1.*z1 = 0, x2.*z2 = 0,
%                     where z = z1 - z2 and z1, z2 > 0.
%                     Typically the same as Featol.
% options.StepTol     (between 0 and 1): Controls how close each an
%                     x or z may be to reaching a bound at each step.
%                     For safety, should not be bigger than 0.99 (say)
%                     for nonlinear problems.
% options.StepSame    1 (true) if stepx and stepz should be the same
%                     (gives true Newton method for nonlinear problems);
%                     0 (false) if stepx and stepz may be different
%                     (gives better performance for linear programs).
% options.x0min       Min distance between x0 and bl or bu  AFTER SCALING.
%                     1.0 is about right for cold starts.
%                     0.1 or 0.01 may be ok if x0, z0 are "good".
% options.z0min       Min distance between abs(z0) and zero AFTER SCALING,
%                     when z0 is used to initialize z1 > 0 and z2 > 0.
%                     Typically the same as x0min.
% options.mu0         Initial mu (ABSOLUTE VALUE) for solving scaled problem.
% options.Method      Specifies how each search direction (dx,dy,dz1,dz2)
%                     should be computed.  Several methods exist for
%                     experimental purposes.  If A has fewer rows than columns
%                     (m < n) we usually solve for dy first (most of the work)
%                     and then get the other components cheaply.
%
%  Method  Solve for  Using                                         Implemented?
%     1       dy      Sparse Cholesky on (A D^2 A' + D2^2 I).           Yes
%     2       dy      Sparse QR on corresponding least-squares problem  Yes
%     3       dy      LSQR on least-squares problem                     Yes
%
%    11       dx      Sparse Cholesky on (D A'A D  + D2^2 I).           No
%    12       dx      Sparse QR on corresponding least-squares problem  No
%    13       dx      LSQR on least-squares problem                     No
%
%    21    dx,dy      Sparse LU on 2x2 KKT-type system                  No
%    23    dx,dy      SYMMLQ on same system                             No
%
%    31    dx,dy,dz1  Sparse LU on 3x3 system (only for problems with   No
%                     with vanilla bounds: 0 < x < inf).
%                     This is a HUGE system, but it is relatively
%                     well-conditioned compared to the other systems.
%
%    41 dx,dy,dz1,dz2 Sparse LU on 4x4 system with general bounds.      Yes
%                     This is an even HUGER system.
%
%    If A is an explicit sparse matrix, all methods are applicable.
%     1 is usually best (e.g. for LPs).
%     2 may be more reliable; it's there for checking.
%     3 is sometimes more efficient (e.g. for entropy problems).
%       Diagonal preconditioning is possible.
%
%    If A is an operator,
%     3 must be used.  Diagonal preconditioning is not possible.
%
%    Notes:
%    Method =  1      On ENTROPY.big, symamd can find an ordering
%                     but chol never finishes.
%    Method = 41      On ENTROPY.big, symamd never finishes.
%
% The following options control LSQR when Method = 3:
%
% options.LSQRMaxIter * min(m,n) is the maximum LSQR (CG) iterations.
% options.LSQRatol1   is the starting value of the LSQR accuracy
%                     tolerance "atol" (if LSmethod = 3).
%                     1e-3 or 1e-4 sometimes works.
%                     1e-8 may be needed for LPs.
%                     In general, if max(Pinf,Dinf,Cinf) doesn't decrease
%                     every iteration, set atol1 to a smaller value.
% options.LSQRatol2   is the smallest value atol is reduced to.
% options.LSQRconlim  shuts LSQR down early if its matrix is ill-conditioned.
%
% options.wait = 0    means pdco should proceed to solve the problem.
%              = 1    means pdco should pause to allow interactive resetting
%                     of some of the parameters.


% pdcoSet.m is derived from optimset.m (Revision 1.14, 1998/08/17)
% in the Optimization Toolbox of The MathWorks, Inc.
%
% 28 Sep 2000: First version of pdscoSet.m.
% 30 Sep 2002: First version of pdcoSet.m derived from pdscoSet.m.
% 04 Oct 2002: StepSame introduced.
% 10 Aug 2003: mu0 is now an absolute value -- the initial mu.
% 19 Sep 2003: Method = 41 implemented.
% 22 Sep 2003: LSproblem and LSmethod replaced by Method.
%
%
% Michael Saunders, SOL, Stanford University.


if (nargin == 0)        % Set default options.
    defoptions.MaxIter      =    30;
    defoptions.FeaTol       =  1e-6;
    defoptions.OptTol       =  1e-6;
    defoptions.StepTol      =  0.99;
    defoptions.StepSame     =     1;  % 1 for stepx == stepz (NLPs)
    defoptions.x0min        =   1.0;  % 1.0 | 0.1 for cold | warm starts?
    defoptions.z0min        =   1.0;  % 1.0 | 0.1 for cold | warm starts?
    defoptions.mu0          =  1e-1;  % < 1.0 better than 1.0?
    defoptions.Method       =     3;  % 3 = computed dy using LSQR
    defoptions.LSQRMaxIter  =  10.0;
    defoptions.LSQRatol1    = 1e-08;
    defoptions.LSQRatol2    = 1e-15;  %
    defoptions.LSQRconlim   = 1e+12;  % Somewhere between e+8 and e+16
    defoptions.wait         =     0;
    defoptions.NOTE         = 'LSQRMaxIter is scaled by the matrix dimension';
    
    if (nargout == 0)    % Display options.
        disp('pdco default options:')
        disp( defoptions )
    else
        options = defoptions;
    end
    return;
end

Names = ...
    [
    'MaxIter    '
    'FeaTol     '
    'OptTol     '
    'StepTol    '
    'StepSame   '
    'x0min      '
    'z0min      '
    'mu0        '
    'Method     '
    'LSQRMaxIter'
    'LSQRatol1  '
    'LSQRatol2  '
    'LSQRconlim '
    'wait       '
    'NOTE       '
    ];
m     = size (Names,1);
names = lower(Names);

% The remaining clever stuff is from optimset.m.
% We should obtain permission from the MathWorks.

% Combine all leading options structures o1, o2, ... in pdcoSet(o1,o2,...).
options = [];
for j = 1:m
    eval(['options.' Names(j,:) '= [];']);
end
i = 1;
while i <= nargin
    arg = varargin{i};
    if isstr(arg)                         % arg is an option name
        break;
    end
    if ~isempty(arg)                      % [] is a valid options argument
        if ~isa(arg,'struct')
            error(sprintf(['Expected argument %d to be a ' ...
                'string parameter name ' ...
                'or an options structure\ncreated with pdcoSet.'], i));
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
                eval(['val = arg.' Names(j,:) ';']);
            else
                val = [];
            end
            if ~isempty(val)
                eval(['options.' Names(j,:) '= val;']);
            end
        end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end

expectval = 0;                          % start expecting a name, not a value

while i <= nargin
    arg = varargin{i};
    
    if ~expectval
        if ~isstr(arg)
            error(sprintf(['Expected argument %d to be a ' ...
                'string parameter name.'], i));
        end
        
        lowArg = lower(arg);
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error(sprintf('Unrecognized parameter name ''%s''.', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches
            % (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
                msg = [msg '(' deblank(Names(j(1),:))];
                for k = j(2:length(j))'
                    msg = [msg ', ' deblank(Names(k,:))];
                end
                msg = sprintf('%s).', msg);
                error(msg);
            end
        end
    else
        eval(['options.' Names(j,:) '= arg;']);
    end
    
    expectval = ~expectval;
    i = i + 1;
end

if expectval
    error(sprintf('Expected value for parameter ''%s''.', arg));
end

%
% Copyright (c) 2006. Michael Saunders
%

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
