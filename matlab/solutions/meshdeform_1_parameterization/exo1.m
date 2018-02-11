q = floor(p/4);
t = (0:(q-1))/q;
t1 = (0:(p-3*q-1))/(p-3*q);
Z = [[t, t*0+1, 1-t, t1*0]; ...
     [t*0, t, t*0+1, 1-t1]];
clf; hold('on');
hh = plot(Z(1,[1:p 1]),Z(2,[1:p 1]), '.-');
if using_matlab()
    set_linewidth(hh,3);
end
hh = plot(Z(1,[1:p 1]),Z(2,[1:p 1]), 'r.');
if using_matlab()
    set_markersize(hh,20);
end
axis off; axis square; 
