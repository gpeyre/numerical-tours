clf;
for i=1:4
    q = qlist(i);
    Dend = Dland( pend(1), pend(2), :);
    H = max(abs(Dland(:,:,1:q)-repmat(Dend(1:q), [n n 1])), [], 3);
    %
    options.end_points = pend;
    options.heuristic = H;
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
    h = plot(landmarks(1,1:q), landmarks(2,1:q), 'y.'); set(h, 'MarkerSize', 15);
    axis ij;
end
