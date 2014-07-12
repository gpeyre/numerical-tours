x = [.5;.5];
niter = 20;
E = [];
D = [];
X = [];
for i=1:niter
    X(:,i) = x;
    E(end+1) = f(x);
    D(end+1) = norm(x);
    x = x - tau*Gradf(x);
end
clf;
h = plot(log10(E));
set(h, 'LineWidth', 2);
axis tight;
title('log_{10}(x^{(k)})');
