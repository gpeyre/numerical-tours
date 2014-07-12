n = 8;
clf;
for k=1:4
    h = hdd(k);
    f = zeros(n,1); f(n/2+1) = 1;
    for i=1:5
        f = subdivide(f,h);
    end
    subplot(4,1,k);
    plot(linspace(-n/2,n/2,length(f)), f); 
    axis([-n/2 n/2 -.15 1.03]);
end
