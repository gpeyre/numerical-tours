function [x,Flist,nbr_linesearch,flist] = NewtonIP(Q,s,C,b,t,x, epsilon)

% solves
%       min_{x in R^d} 1/2*<Q*x,x> + <x,s> - 1/t*sum(log(b-C*x))



dotp = @(x,y)sum(x(:).*y(:));
mylog = @(x)log(x).*(x>0) - (x<=0).*1e20;
f = @(x)1/2*dotp(Q*x,x) + dotp(x,s);
% function penalized with entropy barrier
F = @(x,t)f(x) - 1/t * sum(mylog(b-C*x));
% gradient and hessian
nablaF   = @(x,t)Q*x + s + 1/t*C'*(1./(b-C*x));
hessianF = @(x,t)Q + 1/t*C' * diag(1./(b-C*x).^2) * C;


%%
% Check with finite difference.

if 0
tau = 1e-10; dx = randn(size(x)); 
relerr = @(u,v)norm(u-v)/norm(u);
err = relerr( (F(x+tau*dx,t)-F(x,t))/tau , dotp(nablaF(x,t),dx) );
err = relerr(   dotp(nablaF(x+tau*dx,t)-nablaF(x,t),dx)/tau, ...
                dotp(hessianF(x,t)*dx,dx) );
err = relerr(   (nablaF(x+tau*dx,t)-nablaF(x,t))/tau, ...
                hessianF(x,t)*dx );
end

niter_max_newton = 200; % iteration maximum of newton
niter_line_search = 2000; 
alpha = .3; % error reduction of line search
beta = .7; % step reduction of line search


Flist = [];
flist = [];
nbr_linesearch = [];
for it=1:niter_max_newton
    Flist(it) = F(x,t);
    flist(it) = f(x);
    nabla = nablaF(x,t);
    Delta = -pinv( hessianF(x,t) ) * nabla;
    % check exit
    incr = -dotp( Delta, nabla );
    if incr/2<=epsilon
        break;
    end
    rho = 1; 
    
    if 0
        % just to display the line search
        rholist = linspace(0,1,10000);
        Flist = [];
        for i=1:length(rholist)
            Flist(i) = F(x+rholist(i)*Delta,t);
        end
        Flist(Flist>1e10) = Inf;
        clf; hold on;
        plot(rholist, Flist, 'b');
        plot(rholist, F(x,t) + alpha*rholist*dotp(nabla,Delta), 'r');
    end
    
    % line search  
    for jt=1:niter_line_search
        if F(x+rho*Delta,t)<F(x,t) + alpha*rho*dotp(nabla,Delta)
            break; 
        end
        rho = beta*rho;
    end
    nbr_linesearch(it) = jt;
    if jt==niter_line_search
        warning('Maximum line search reached');
    end
    % move
    x = x + rho*Delta;    
end
if it==niter_max_newton
    warning('Maximum Newton iteration reached');
end

end