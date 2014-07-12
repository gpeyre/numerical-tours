niter = 300;
x = y;
E = [];
for i=1:niter
    E(end+1) = f(y,x,epsilon);
    x = x - tau*Gradf(y,x,epsilon);
end
clf;
h = plot(E);
set(h, 'LineWidth', 2);
axis tight;
