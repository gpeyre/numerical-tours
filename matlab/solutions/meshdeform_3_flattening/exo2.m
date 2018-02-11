niter = 150;
stress = [];
vertexS = vertexF;
ndisp = [1 5 10 min(niter,100) Inf];
k = 1;
clf;
for i=1:niter
    if ndisp(k)==i
        subplot(2,2,k);
        plot_mesh(vertexS,faces,options);
        k = k+1;
    end
    % Compute the distance matrix.
    D1 = repmat(sum(vertexS.^2,1),n,1);
    D1 = sqrt(D1 + D1' - 2*vertexS'*vertexS);
    % Compute the scaling matrix.
    B = -D./max(D1,1e-10);
    B = B - diag(sum(B));
    % update
    vertexS = (B*vertexS')' / n;
    % Xstress = Xstress-repmat(mean(Xstress,2), [1 n]);
    % record stress
    stress(end+1) = sqrt( sum( abs(D(:)-D1(:)).^2 ) / n^2 );
end
