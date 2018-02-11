p = round(n/3);
faces1 = faces;
vertex1 = vertex;
for i=1:n-p
    edges = compute_edges(faces1);
    nedges = size(edges,2);
    k = floor(rand*(nedges-1))+1;
    e = edges(:,k);
    vertex1(:,e(1)) = mean( vertex1(:,e),2 );
    vertex1(:,e(2)) = Inf;
    faces1(faces1==e(2)) = e(1);
    a = sum( diff(sort(faces1))==0 );
    faces1(:,a>0) = [];
end
% display
clf;
plot_mesh(vertex1,faces1,options);
shading faceted;
