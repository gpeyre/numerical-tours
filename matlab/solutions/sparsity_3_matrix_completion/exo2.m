k = 30; % number of trial per sparsity
rmin = round(.5*P/(n*log(n)));
rmax = round(2*P/(n*log(n)));
rlist = rmin:2:rmax;
niter = 200;
R = zeros(length(rlist),1);
tol = .02;
for ir=1:length(rlist)
    % progressbar(ir,length(rlist));
    r = rlist(ir);
    for j=1:k
        x0 = randn(n,r)*randn(r,n);
        y = Phi(x0);
        ProxF = @(x,gamma)x + PhiS(y-Phi(x)); 
        rProxF = @(x,gamma)2*ProxF(x,gamma)-x;
        %
        x = x0; tx = x;
        for i=1:niter
            tx = (1-mu/2)*tx + mu/2*rProxG( rProxF(tx,gamma),gamma );
            x = ProxF(tx,gamma);
        end
        R(ir) = R(ir) + (norm(x-x0)/norm(x0)<tol);
    end    
end
clf;
h = plot(rlist, R/k); 
axis([min(rlist) max(rlist) -.05 1.05]);
xlabel('r');
set(h, 'LineWidth', 2);
