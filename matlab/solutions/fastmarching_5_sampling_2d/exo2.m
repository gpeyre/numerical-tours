nmax = 200;
ndisp = [20 50 100 200];
vertex = [1;1];
[D,Z,Q] = perform_fast_marching(1./W, vertex);
clf; u = 1;
vertex_svg = {};
faces_svg = {};
for k=2:nmax
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
        u = u+1;  
        % compute the Delaunay triangulation
        [D1,Z,Q] = perform_fast_marching(1./W, vertex);
        vertex_svg{end+1} = vertex;
        faces_svg{end+1} = compute_voronoi_triangulation(Q,vertex);
    end
end
