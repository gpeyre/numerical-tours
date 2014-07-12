wlist = [.6 .8 .9 .95];
clf;
for i=1:4
    weight = wlist(i);
    options.end_points = pend;
    options.heuristic = weight*H;
    options.nb_iter_max = Inf;
    options.constraint_map = Inf+zeros(n);
    [D,S] = perform_fast_marching(1./W, pstart, options);
    %
    I = find(S<0);
    U = cat(3,M,M,M);
    U(I) = 1; U([I+n^2, I+2*n^2]) = U([I+n^2, I+2*n^2])*.3;
    subplot(2,2,i); 
    hold on;
    imageplot(U);
    h = plot(p(2,:),p(1,:), '.k'); set(h, 'LineWidth', 2);
    h = plot(pstart(2),pstart(1), '.g'); set(h, 'MarkerSize', 25);
    h = plot(pend(2),pend(1), '.b'); set(h, 'MarkerSize', 25);
    axis ij;
end
