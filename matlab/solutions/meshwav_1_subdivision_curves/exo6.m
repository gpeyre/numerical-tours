w = 1/16;
h = h4pt(w);
f = f0;
clf;
for j=0:Jmax
    f = subdivide(f,h);
    subplot(2,2,j+1);
    hold on;
    myplot(f, 'k.-');
    myplot(f0, 'r.--');    
    myaxis(.13);
    hold off;
end
