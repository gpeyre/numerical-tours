function [xtv,err,tv,lalist,Err] = perform_tv_denoising(x,options)

% perform_tv_denoising - denoising with TV minimization
%
%   [xtv,err,tv,lalist,Err] = perform_tv_denoising(x,options);
%
%   Solve the lagragian TV minimization (Rudin-Osher-Fatemi)
%       xtv = argmin_u 1/2 * |x-u|^2 + lambda*TV(u)
%   where TV is the discrete TV norm
%       1D:  TV(u) = sum_i |u(i)-u(i-1)|
%       2D:  TV(u) = sum_i sqrt( ux(i)^2 + uy(i) )
%   where ux,uy are partial derivative (finite difference along x/y).
%
%   It implements Chambolle's iterative fixed point algorithm
%       A. Chambolle: An Algorithm for Total Variation Minimization and Applications.
%       Journal of Mathematical Imaging and Vision 20 (1-2): 89-97, January - March, 2004
%
%   Actually, the original Chambolle's impementation is done when
%   options.method='chambolle', and a projected gradient iteration is
%   applied if options.method='gradient'. The projected gradient iteration
%   seems to be slightly faster.
%
%   Works for 1D signal and 2D images.
%
%   lambda controls how much denoising you want.
%   You can track the denoising level in various way:
%       * simply set options.lambda
%       * set values options.lambda_max and options.lambda_min
%       and the algorithm will linearely decrease the value of 
%       lambda between these values.
%       If you set (in an oracle fashion) options.x0, then 
%       the algorihtm will return the xtv that is the closest to x0.
%       * set a target error tolerance options.etgt and the algorithm
%       will return the xtv such that norm(xtv-x,'fro')=etgt
%       * set a target TV norm options.tvtgt and the algorithm will
%       return xtv such that TV(xtv)=tvtgt.
%
%   To perform TV inpainting, you can set options.mask and options.xtgt
%   where mask(i)=1 is a missing pixel and mask(i)=0 impose value xtgt(i)
%   for the solution.
%
%   Copyright (c) 2008 Gabriel Peyre


options.null = 0;

nbdims = nb_dims(x);


la = getoptions(options, 'lambda', []);
lamax = .4;
lamin = eps;
if not(isempty(la))
    lamax = la; lamin = la;
end
lamax = getoptions(options, 'lambda_max', lamax);
lamin = getoptions(options, 'lambda_min', lamin);
verb = getoptions(options, 'verb', 1);
lamin = max(lamin,eps);
method = getoptions(options, 'method', 'chambolle');

% possible use for inpainting
mask = getoptions(options, 'mask', []);
xtgt = getoptions(options, 'xtgt', []);

niter = getoptions(options, 'niter', 3000);
% max iteration for inner loop
niter_inner = getoptions(options, 'niter_inner', 30);
% tolerance for early stop in inner loop
tol = getoptions(options, 'tol', 1e-5);

% use for error tracking
etgt = getoptions(options, 'etgt', []);
% use for TV tracking
tvtgt = getoptions(options, 'tvtgt', []);
if not(isempty(etgt)) && not(isempty(tvtgt))
    warning('Cannot track bot error and TV.');
    tvtgt = [];
end

x0 = getoptions(options, 'x0', []);
tau = getoptions(options, 'tau', 2/8);
if tau>2/8
    warning('tau is too large for convergence');
    tau = 2/8;
end
if nbdims==1
    xi = zeros(size(x))
else
    xi = zeros([size(x) 2]);
end
x1 = x;
xi = getoptions(options,'xi', xi );
display = getoptions(options,'display', 1);
TVoperator = getoptions(options, 'TVoperator', []);
if not(isempty(TVoperator))
    y = getoptions(options, 'TVrhs', [], 1); % mendatory
end

if isempty(TVoperator) && isempty(mask)
    niter_inner = 1;
end

lalist = [];
if isempty(etgt)
    lalist = linspace(lamax,lamin,niter);
end

if (not(isempty(etgt)) || not(isempty(tvtgt))) && isempty(la)
    la = .1;
end

nrefresh = getoptions(options, 'nrefresh', 3);  % rate of refresh for lambda tracking
ndisp = max(round( niter/100 ),2);

if display
    clf;
end
err = []; tv = []; Err = [];
for i=1:niter
    
    if verb
        progressbar(i,niter);
    end
    
    if isempty(etgt) && isempty(tvtgt)
        % fixed or decaying lambda
        la = lalist(i);
    end
    
    if not(isempty(xtgt)) && not(isempty(mask)) % && i>niter/4
        % inpainting mode
        x = x - la*div( xi, options );
        x(mask==0) = xtgt(mask==0);
    end
    if not(isempty(TVoperator))
        % perform projection on constraintes
        x = x - la*div( xi, options );
        r = y - feval(TVoperator, x, +1, options);
        err(i) = norm(r, 'fro');
        x = x + feval(TVoperator, r, -2, options);
        x = real(x);
    end
    
    %%% INNER LOOP
    for iinner = 1:niter_inner
        gdv = grad( div(xi,options) - x/la, options );
        if strcmp(method, 'chambolle')
            %% CHAMBOLLE STEP
            d = sqrt(sum(gdv.^2,3));
            d = repmat( d, [1 1 nbdims] );
            xi = ( xi + tau*gdv ) ./ ( 1+tau*d );
        else
            %% GRADIENT STEP
            xi = xi + tau*gdv;
            d = sqrt(sum(xi.^2,3)); d = max(d,1e-10); 
            d = repmat( d, [1 1 nbdims] );
            xi = xi .* min(d,1)./d;
        end
        % reconstruct
        x2 = x - la*div( xi, options );
        % early stop ?
        relerr = norm(x2-x1,'fro')/sqrt(prod(size(x1)));
        if relerr<tol
            x1 = x2; break;
        end
        x1 = x2; % new solution
    end    
    
    if nargout>=2
        if isempty(TVoperator)
            err(i) = norm( x-x1, 'fro' );
        end
    end
    if nargout>=3
        tv(i) = compute_total_variation(x1,options);
    end

    if mod(i,nrefresh)==0 && (not(isempty(etgt)) || not(isempty(tvtgt)))
        % update value of lambda
        if not(isempty(etgt))
            mu = etgt / norm( x-x1, 'fro' );
        end
        if not(isempty(tvtgt))
            mu = compute_total_variation(x1) / tvtgt;
        end
        mu = max(min(mu,2),1/2);
        la = la*mu;
        lalist(i) = la;
    end

    if mod(i,ndisp)==1 && display==1
        if nbdims==1
            plot(x1); axis tight; % axis([1 n -eta 1+eta]);
        else
            imageplot(x1); 
        end
        drawnow;
    end
     
    if not(isempty(x0))
        % oracle error
        Err(end+1) = norm(x0-x1, 'fro');
    end
end
xtv = x1;

function d = nb_dims(x)

% nb_dims - debugged version of ndims.
%
%   d = nb_dims(x);
%
%   Copyright (c) 2004 Gabriel Peyré

if isempty(x)
    d = 0;
    return;
end

d = ndims(x);

if d==2 && (size(x,1)==1 || size(x,2)==1)
    d = 1;
end