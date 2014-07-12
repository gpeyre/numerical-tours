alpha_list = [.4 .7 1 1.3];
m = 600;
clf;
for ialpha=1:length(alpha_list)
    alpha = alpha_list(ialpha);
    % metric
    W = sum(G.^2,3);
    W = perform_blurring(W, 6);
    W = (W+epsilon).^alpha;
    % Farthest point
    vertex = [1;1];
    [D,Z,Q] = perform_fast_marching(1./W, vertex);
    for k=2:m
        options.constraint_map = D;
        [D1,Z,Q] = perform_fast_marching(1./W, vertex(:,end), options);
        D = min(D,D1);
        [tmp,i] = max(D(:));
        [x,y] = ind2sub([n n],i);
        vertex(:,end+1) = [x;y];
    end
    % compute the Delaunay triangulation
    [D,Z,Q] = perform_fast_marching(1./W, vertex);
    faces = compute_voronoi_triangulation(Q,vertex);
    % compute approximation
    vgeod = compute_orthoproj_triangulation(vertex, faces, M);
    Mgeod = compute_triangulation_interpolation(faces,vertex,vgeod, n);
    % display
    imageplot(clamp(Mgeod), ['\alpha= ' num2str(alpha) ', SNR=' num2str(snr(Mgeod,M),3) 'dB'], 2,2,ialpha);
end
