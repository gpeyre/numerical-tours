face = face0;
vertex = vertex0;
clf;
for i=1:4
    edge = compute_edges(face);
    n = size(vertex,2);
    ne = size(edge,2);
    % Compute the number of the three edges associated to each face.
    A = sparse([edge(1,:);edge(2,:)],[edge(2,:);edge(1,:)],[n+(1:ne);n+(1:ne)],n,n);
    v12 = full( A( face(1,:) + (face(2,:)-1)*n ) );
    v23 = full( A( face(2,:) + (face(3,:)-1)*n ) );
    v31 = full( A( face(3,:) + (face(1,:)-1)*n ) );
    % Compute the new faces, each old face generates 4 faces.
    face = [  cat(1,face(1,:),v12,v31),...
        cat(1,face(2,:),v23,v12),...
        cat(1,face(3,:),v31,v23),...
        cat(1,v12,v23,v31)   ];
    % Add new vertices at the edges center.
    vertex = [vertex, (vertex(:,edge(1,:))+vertex(:,edge(2,:)))/2 ];
    % Project the points on the sphere.
    d = sqrt( sum(vertex.^2,1) );
    vertex = vertex ./ repmat( d, [size(vertex,1) 1]);
    % display
    subplot(2,2,i);
    plot_mesh(vertex,face);
    shading('faceted'); lighting('flat'); view(3); axis('tight');
end
