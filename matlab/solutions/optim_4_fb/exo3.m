gamma = 1/L;
x = y;
z = y;
t = 1;
Efista = [];
for i=1:niter
    xold = x;
    told = t;
    x = proxg( z - gamma*gradf(z), gamma );
    t = ( 1 + sqrt(1+4*t^2) )/2;
    z = x + (told - 1)/t*( x-xold );
    % record the energy degay.
    Efista(i) = 1/2*norm(Phi(x) - y)^2 + lambda*norm(x, 1);     
end
% display energy decay
Emin = min([min(E) min(Efista) min(Erelax(:))]);
clf; hold on;
hh = plot(log10( E(1:end/4)-Emin ), 'b');
set(hh, 'LineWidth', 2);
hh = plot(log10( Erelax(1:end/4,end)-Emin ), 'g');
set(hh, 'LineWidth', 2);
hh = plot(log10( Efista(1:end/4)-Emin ), 'r');
set(hh, 'LineWidth', 2);
axis('tight');
set_label('i', 'log_{10}( E-E* )');
legend('FB', ['FB-relax, \mu=' num2str(mu_list(end))], 'FISTA');
