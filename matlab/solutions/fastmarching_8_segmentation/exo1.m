v = sort(D(:));
Tlist = v( round([.05 .1 .15 .25]*n^2) );
clf;
for i=1:4
    T = Tlist(i);
    A = repmat(M,[1 1 3]);
    I = find(D<T);
    A(I) = 1;
    A([I+n^2;I+2*n^2]) = 0;
    subplot(2,2,i);
    hold on;
    imageplot(A);
    [c,h] = contour(D<T, [.5 .5], 'b');
    set(h, 'LineWidth', 2.5); 
    axis('ij');
end
