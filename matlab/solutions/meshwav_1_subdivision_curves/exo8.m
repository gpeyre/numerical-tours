hh = [];
hh{end+1}=[1 3 3 1]/4;
hh{end+1} = [1 4 6 4 1]/8;
hh{end+1} = [-1, 0, 9, 1, 9, 0, -1]/16; 
hh{end}((end+1)/2)=1;
lgd = {'Quadratic', 'Cubic', 'Interpolating'};
col = {'r' 'g' 'b'};    
clf; hold on;
for k=1:length(hh)
    h = hh{k};
    f = f0;
    for j=0:7
        f = subdivide(f, h);
    end
    myplot(f, col{k});
end
myplot(f0, 'k.--');
axis tight; axis off;
legend(lgd);
