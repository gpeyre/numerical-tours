mu_list = [-.5 0 .5 .8 .9 .95 .99];
Erelax = [];
for imu = 1:length(mu_list)
    mu = mu_list(imu);
    lgd{imu} = ['\mu=' num2str(mu)];
    x = y; z = y;
    for i=1:niter
        xold = x;
        told = t;
        x = proxg( z - gamma*gradf(z), gamma );
        z = x + mu*( x-xold );
        % record the energy degay.
        Erelax(i,imu) = 1/2*norm(Phi(x) - y)^2 + lambda*norm(x, 1);
    end
end
% display
Emin = min(Erelax(:));
clf; hold on;
hh = plot(log10( Erelax(1:end/4,:)-Emin ) );
set(hh, 'LineWidth', 2);
axis('tight');
set_label('i', 'log_{10}( E-E* )');
legend(lgd);
