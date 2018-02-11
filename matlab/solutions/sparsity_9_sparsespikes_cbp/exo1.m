A = @(u)GammaS(Gamma(u));
%
u = randn([N 2]);
u = u/norm(u(:));
%
e = [];
for i=1:15
    v = A(u);
    e(end+1) = sum(u(:).*v(:));
    u = v/norm(v(:));
end
L = e(end);
clf; plot(e, 'LineWidth', 2); axis tight;
