function [x,err,it] = perform_conjugate_gradient(A,y,options)


% perform_conjugate_gradient - perform (bi)-conjugate gradient
%
%     [x,err,k] = perform_conjugate_gradient(A,y,options);
%
%   Solves for A*x=y.
%   Works both for vector x,y and matrix x,y (parralel soving of many
%   linear equations).
%
%   Important: the algorithm assumes that the matrix is symmetric definite
%   positive. If it is not the case, then you must set options.is_sdp=0,
%   and the algorithm will use a bi-conjugate gradient descent (x2 slower).
%
%   A can be a matrix or a callback function y=A(x,options).
%   In this case, you have to set options.ncols as the number of columns of
%   A (if it is a callback).
%   In the non-symmetric case, then options.transpose=1 for the computation
%   of A*x and options.transpose=-1 for the computation of A'*x.
%
%   err monitors the decreasing of |A*x-y|.
%   k is the total number of iterations.
%
%   You can set:
%       options.x is an initial guess
%       options.epsilon is maximum error
%       options.niter_max is maximum number of iterations
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
niter = getoptions(options, 'niter_max', 100);
epsilon = getoptions(options, 'epsilon', 1e-5);
is_sdp = getoptions(options, 'is_sdp', 1);
x = getoptions(options, 'x', []);

use_callback = 0;
if not(isnumeric(A))
    use_callback = 1;
end
if not(isempty(x))
    x = options.x;
else
    if use_callback==0
        x = zeros(size(A,2),1);
    else
        if isfield(options, 'ncols')
            x = zeros(options.ncols, 1);
        else
            error('You have to specify options.ncols');
        end
    end
end


if norm(y(:), 'fro') == 0
    normb = epsilon;
else
    normb = epsilon; % * sum(y(:).^2);
end

if use_callback==0
    r  = y - A*x;
else
    options.transpose = 1;
    r  = y - feval(A,x,options);    
end
p = r;
if is_sdp==0
    rr = r;
    pp = p;
end
r0 = sum(r(:).^2);


err = [sum(r0)];
for it=1:niter
    
    % auxiliary vector
    if use_callback==0
        w  = A*p;
    else
        options.transpose = 1;
        % options.transpose = -1;
        w  = feval(A,p,options);
    end
    
    if is_sdp==1
        d = sum(p(:) .* w(:));
    else
        d = sum(pp(:) .* w(:));
    end
    I = find(abs(d)<eps); d(I) = 1;
    alpha = repmat( r0 ./ d, [size(x,1) 1] );              % found optimal alpha in line search
    x = x + alpha.*p;                       % new guess
    r = r - alpha.*w;                       % the residual is in fact r=b-A*x
    if is_sdp==0
        if use_callback==0
            w  = A'*pp;
        else
            options.transpose = -1;
            w  = feval(A,pp,options);
        end
        rr = rr-alpha.*w;
    end
    
    rprev = r0;                             % save norm of the old residual
    if is_sdp==0
        r0 = sum(r.*rr);
    else
        r0 = sum(r.^2);                         % compute norm of new residual
    end
    
    err(end+1) = sqrt( sum(r.^2) );
    if err(end)<normb
        return;
    end

    % search direction
    beta = r0./rprev;
    p = r + repmat(beta, [size(x,1) 1]).*p;
    if is_sdp==0
        pp = rr + repmat(beta, [size(x,1) 1]).*pp;
    end


end

