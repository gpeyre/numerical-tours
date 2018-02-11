niter = 4000;
x = y;
E = [];
for i=1:niter
    x = proxg( x - gamma*gradf(x), gamma );
    % record the energy degay.
    E(i) = 1/2*norm(Phi(x) - y)^2 + lambda*norm(x, 1);     
end
% display energy decay
clf;
hh = plot(log10( (E(1:end/4)-E(end))/(E(1)-E(end)) ));
set(hh, 'LineWidth', 2);
axis('tight');
set_label('i', 'log_{10}( E-E* )');
