m = 200;
% initialize
landmarks = [100];
[D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks);
clf;
k = 1; displist = [100 m];
for i=2:m
    % select
    [tmp,landmarks(end+1)] = max(D);
    % update
    options.constraint_map = D;
    [D1,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks,options);
    D = min(D,D1);
    if i==displist(k)
        [D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks);
        % compute the mesh
        V = Q(faces); V = sort(V,1);
        V = unique(V', 'rows')';
        d = 1 + (V(1,:)~=V(2,:)) + (V(2,:)~=V(3,:));
        %
        I = find(d==3); I = sort(I);
        z = zeros(n,1);
        z(landmarks) = (1:length(landmarks))';
        facesV = z(V(:,I));
        vertexV = vertex(:,landmarks);
        % Re-orient the faces so that they point outward of the mesh.
        options.method = 'slow';
        options.verb = 0;
        facesV = perform_faces_reorientation(vertexV,facesV, options);
        % display
        subplot(1,2,k);
        options.face_vertex_color = [];
        plot_mesh(vertexV,facesV, options);
        shading faceted;
        %
        k = k+1;
    end
end
