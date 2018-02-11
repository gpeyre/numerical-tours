nlist = round(linspace(100,4000, 10));
options.niter = 200;
err = [];
for i=1:length(nlist)
    n = nlist(i);    
    warning off; rand('state', 0); warning on;
    vertex = [2*rand(2,n)-1; zeros(1,n)];
    vertex(:,1) = 0;
    faces = compute_delaunay(vertex);
    options.U = sqrt(sum(vertex.^2))';
    [U,e] = perform_geodesic_iterative(vertex, faces, W, I, options);
    err(end+1) = norm(U-options.U)/sqrt(n);
end
clf;
h = plot(nlist, err, '.-'); 
set(h, 'LineWidth', 2);
axis tight;
