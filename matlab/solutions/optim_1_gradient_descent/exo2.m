taulist = [.3 1 1.7]/eta;
xinit = [[.7;.7],[-.7;.5],[-.7;-.6]];
collist = {'k' 'g' 'r'};
%
clf; hold on;
imagesc(t,t,F); colormap jet(256);
contour(t,t,F, 20, 'k');
for k=1:length(taulist);
    tau = taulist(k);
    %
    x = xinit(:,k);
    niter = 100;
    X = [];
    for i=1:niter
        X(:,i) = x;
        x = x - tau*Gradf(x);
    end
    %
    h = plot(X(1,:), X(2,:), [collist{k} '.-']);
    set(h, 'LineWidth', 2);
    set(h, 'MarkerSize', 15);
    axis off; axis equal;
end
