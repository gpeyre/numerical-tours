p = round(n/3);
faces1 = faces;
vertex1 = vertex;
for i=1:n-p
    edges = compute_edges(faces1);
    D = vertex(:,edges(1,:)) - vertex(:,edges(2,:));
    D = sum(D.^2,1);
    [tmp,k] = min(D);
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
