d = (D(pend(1),pend(2)) + H(pstart(1),pstart(2)))/2;
Tlist = [1.01 1.05 1.1 1.2] * d;
clf;
for i=1:4
    t = Tlist(i);
    U = cat(3,M,M,M);
    I = find( D+H<=t );
    U(I) = 1; U([I+n^2, I+2*n^2]) = U([I+n^2, I+2*n^2])*.3;
    subplot(2,2,i);
    hold on;
    imageplot(U);
    h = plot(p(2,:),p(1,:), '.k'); set(h, 'LineWidth', 2);
    h = plot(pstart(2),pstart(1), '.g'); set(h, 'MarkerSize', 25);
    h = plot(pend(2),pend(1), '.b'); set(h, 'MarkerSize', 25);
    axis ij;
end
