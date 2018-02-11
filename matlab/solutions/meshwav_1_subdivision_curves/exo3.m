f0 = [0+0i; 1+0i; 1+1i; 0+1i];
f = f0;
clf;
for j=0:Jmax
    f = subdivide(f,h);
    subplot(2,2,j+1);
    hold on;
    myplot(f, 'k.-');
    myplot(f0, 'r.--');
    myaxis(.03);
end
