tf = double( w>0 );
tu = Grad(tf);
E = [];
lun = @(u)sum( sum( sqrt(sum(u.^2,3)) ) );
for i=1:niter
    tf1 = tf; tu1 = tu;
    % rProxG
    [tf2,tu2] = ProxG(tf1, tu1, gamma);
    tf1 = 2*tf2 - tf1; tu1 = 2*tu2 - tu1;
    % rProxF
    [tf2,tu2] = ProxF(tf1, tu1, gamma);
    tf1 = 2*tf2 - tf1; tu1 = 2*tu2 - tu1;
    % Average
    tf = (1-mu/2)*tf + mu/2*tf1;
    tu = (1-mu/2)*tu + mu/2*tu1;
    % Converging sequence
    [f,u] = ProxG(tf,tu,gamma);
    % 
    E(end+1) = sum(w(:).*f(:)) + lambda*lun(u);
end
clf;
plot(log10(E(1:end/2)-min(E)));
axis tight;
title('log_{10}(E(f) - E^*)');
