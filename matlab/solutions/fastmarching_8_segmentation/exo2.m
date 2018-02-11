col = {'r' 'g' 'b' 'c' 'y' 'k' 'r:'};
clf;
hold on;
imageplot(M);
Vlist = unique(Q(:))';
for i=Vlist
    U = zeros(n); U(Q==i)=1;
    [c,h] = contour(U, [.5 .5], col{i});
    set(h, 'LineWidth', 4); 
end
axis ij;
hold off;
