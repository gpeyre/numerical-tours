L = compute_mesh_laplacian(vertex,faces,'combinatorial',options);
[U,S] = eig(full(L)); S = diag(S);
[S,I] = sort(S,'ascend'); U = U(:,I);
vertexF = U(:,2:3)';
vertexF = vertexF - repmat(vertexF(:,icenter), [1 n]);
theta = -pi/2+atan2(vertexF(2,irotate),vertexF(1,irotate));
vertexF = [vertexF(1,:)*cos(theta)+vertexF(2,:)*sin(theta); ...
           -vertexF(1,:)*sin(theta)+vertexF(2,:)*cos(theta)];
clf;
plot_mesh(vertexF,faces);
