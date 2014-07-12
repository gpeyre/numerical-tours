niter = 150;
stress = [];
Xstress = X;
ndisp = [1 5 10 min(niter,100) Inf];
k = 1;
clf;
for i=1:niter
    if ndisp(k)==i
        subplot(2,2,k);
        hold on;
        scatter3(Xstress(1,:),Xstress(2,:),Xstress(3,:),ms,v, 'filled');
        plot_graph(A, Xstress, options);
        colormap jet(256);
        view(v1,v2); axis('equal'); axis('off');
        k = k+1;
    end
    % Compute the distance matrix.
    D1 = repmat(sum(Xstress.^2,1),n,1);
    D1 = sqrt(D1 + D1' - 2*Xstress'*Xstress);
    % Compute the scaling matrix.
    B = -D./max(D1,1e-10);
    B = B - diag(sum(B));
    % update
    Xstress = (B*Xstress')' / n;
    % Xstress = Xstress-repmat(mean(Xstress,2), [1 n]);
    % record stress
    stress(end+1) = sqrt( sum( abs(D(:)-D1(:)).^2 ) / n^2 );
end
