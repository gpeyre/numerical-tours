ndisp = [20 50 100 m];
vertex = [1;1];
[D,Z,Q] = perform_fast_marching(1./W, vertex);
vertex_svg = {};
faces_svg = {};
clf; u = 1;
for k=2:m
    options.constraint_map = D;
    [D1,Z,Q] = perform_fast_marching(1./W, vertex(:,end), options);
    D = min(D,D1);
    [tmp,i] = max(D(:));
    [x,y] = ind2sub([n n],i); 
    vertex(:,end+1) = [x;y];
    if k==ndisp(u)
        subplot(2,2,u);
        hold on;
        imageplot(M, [num2str(k) ' points']);
        plot(vertex(2,:), vertex(1,:), 'r.');
        axis ij;
        u = u+1;  
    end
end
% compute the Delaunay triangulation
[D,Z,Q] = perform_fast_marching(1./W, vertex);
faces = compute_voronoi_triangulation(Q,vertex);
