n = 6;
f = zeros(n,1); f(n/2+1) = 1;
for i=1:5
    f = subdivide(f,h);  
end
clf;
plot(linspace(-1/2,1/2,length(f)), f); axis([-1/2 1/2 -.01 max(f)*1.03]);
