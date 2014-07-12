p = 28;
t0 = (0:1/n:1-1/n)';
t = (0:1/p:1-1/p)';
f0 = interp1(t0,F,t);
f = f0;
Jmax = ceil(log2(n/p));
for j=0:Jmax
    f = subdivide(f, h);
end
clf; hold on;
myplot(F, 'k'); 
myplot(f0, 'k.');
myplot(f, 'r'); 
myaxis(0);
