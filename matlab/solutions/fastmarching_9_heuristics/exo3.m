qlist = [2 5 10 50];
q = max(qlist);
landmarks = floor(rand(2,q)*n)+1;
Dland = zeros(n,n,q);
for i=1:q
    Dland(:,:,i) = perform_fast_marching(1./W, landmarks(:,i));
end
clf;
for i=1:4
    q = qlist(i);
    Dend = Dland( pend(1), pend(2), :);
    H = max(abs(Dland(:,:,1:q)-repmat(Dend(1:q), [n n 1])), [], 3);
    subplot(2,2,i);
    hold on;
    imageplot(H);
    contour(H, 10, 'k', 'LineWidth', 2);
    colormap jet(256);
    h = plot(landmarks(1,1:q), landmarks(2,1:q), 'y.');
    set(h, 'MarkerSize', 15);
    axis ij;
end
