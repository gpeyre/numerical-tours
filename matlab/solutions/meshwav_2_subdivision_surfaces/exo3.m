face = face0;
vertex = vertex0;
for i=1:2
    % Perform the subdivision.
    edge = compute_edges(face);
    n = size(vertex,2);
    ne = size(edge,2);
    % Compute the number of the three edges associated to each face.
    A = sparse([edge(1,:);edge(2,:)],[edge(2,:);edge(1,:)],[n+(1:ne);n+(1:ne)],n,n);
    v12 = full( A( face(1,:) + (face(2,:)-1)*n ) );
    v23 = full( A( face(2,:) + (face(3,:)-1)*n ) );
    v31 = full( A( face(3,:) + (face(1,:)-1)*n ) );
    % Compute the new faces, each old face generates 4 faces.
    face_old = face;
    face = [  cat(1,face(1,:),v12,v31),...
        cat(1,face(2,:),v23,v12),...
        cat(1,face(3,:),v31,v23),...
        cat(1,v12,v23,v31)   ];
    % Compute the vertex and face ring.
    global vring e2f fring facej;
    vring = compute_vertex_ring(face);
    e2f = compute_edge_face_ring(face_old);
    fring = compute_face_ring(face_old);
    facej = face_old;
    % Compute the interpolated position using
    for k=n+1:n+ne
        [e,v,g] = compute_butterfly_neighbors(k, n);
        vertex(:,k) = 1/2*sum(vertex(:,e),2) + 1/8*sum(vertex(:,v),2) - 1/16*sum(vertex(:,g),2);
    end
end
clf;
subplot(1,2,1);
plot_mesh(vertex0,face0,options); 
shading('faceted'); lighting('flat'); axis('tight');
subplot(1,2,2);
plot_mesh(vertex,face,options); 
shading('interp'); lighting('phong'); axis('tight');
