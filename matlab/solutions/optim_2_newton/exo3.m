niter = 12; % 300;
x = y;
E = [];
d = zeros(n);
for i=1:niter
    E(end+1) = f(x);
%    d = Hinv(Gradf(x), A(x), d);
    [d,~] = Hinv(Gradf(x), A(x), d);
    d = flatI(d);
    x = x - d;
end
clf;
h = plot(E);
title('f(x^{l})');
set(h, 'LineWidth', 2);
axis tight;
