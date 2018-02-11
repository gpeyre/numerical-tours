Jmax = 3;
f = f0;
clf;
for j=0:Jmax
    f = subdivide(f,h);
    subplot(2,2,j+1);
    hold on;
    myplot(f, 'k.-');
    myplot(f0, 'r.--');
    myaxis(0);
end
