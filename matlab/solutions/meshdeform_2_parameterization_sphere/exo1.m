vertex1 = vertex;
vertex1 = vertex1 - repmat( mean(vertex1,2), [1 n] );
vertex1 = vertex1 ./ repmat( sqrt(sum(vertex1.^2,1)), [3 1] );
niter = 500;
ndisp = round([0 0.1 0.3 1]*niter); ndisp = max(ndisp,1);
k = 1;
clf;
ninvert = [];
Edir = [];
for i=1:niter
    % smooth 
    vertex1 = vertex1*tW';
    % project
    vertex1 = vertex1 ./ repmat( sqrt(sum(vertex1.^2,1)), [3 1] );
    % record Dirichlet energy
    E = 0;
    for j=1:3
        j1 = mod(j,3)+1;
        % directed edge
        u = vertex1(:,faces(j,:)) - vertex1(:,faces(j1,:));
        % norm squared
        u = sum(u.^2);
        % weights between the vertices
        w = W(faces(j,:) + (faces(j1,:)-1)*n);
        E = sum( w.*u );
    end
    Edir(end+1) = E;
    % record number of inverted triangles
    [normal,normalf] = compute_normal(vertex1,faces);
    C = squeeze(mean(reshape(vertex1(:,faces),[3 3 m]), 2));
    I = sum(C.*normalf);
    ninvert(end+1) = sum(I<0);
    if i==ndisp(k)
        % display
        subplot(2,2,k);
        options.face_vertex_color = double(I(:)>0);
        plot_mesh(vertex1,faces,options);
        colormap gray(256); axis tight;
        shading faceted;
        k = k+1;
    end
end
