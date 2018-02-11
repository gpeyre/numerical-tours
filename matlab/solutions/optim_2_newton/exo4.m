%
tau = 1.8/( 1 + lambda*8/epsilon );
tau = tau*4;
x = y;
E1 = [];
for i=1:niter
    E1(end+1) = f(x);
    x = x - tau*Gradf(x);
end
clf;
h = plot([E;E1]');
title('f(x^{l})');
legend('Newton', 'Gradient');
set(h, 'LineWidth', 2);
axis tight;
