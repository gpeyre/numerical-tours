R = [];
u = zeros(N,2);
niter = 2000;
damp = 1.8;
tau = damp/L;
for i=1:niter
    u = ProxJ( u - tau * gradF(u), lambda*tau );
    R(end+1) = E(u); 
end
sel = 1:niter/4;
plot(sel, log(R(sel)/min(R)-1), '-', 'LineWidth', 2); axis tight;
