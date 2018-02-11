vertex1 = vertex; 
vertex1 = vertex1 - repmat( mean(vertex1,2), [1 n] );
vertex1 = vertex1 ./ repmat( sqrt(sum(vertex1.^2,1)), [3 1] );
niter = 500;
eta = .5;
ndisp = round([0 0.1 0.3 1]*niter); ndisp = max(ndisp,1);
k = 1;
Edir = [];
for it=1:niter
    % compute the center
    A = squeeze(sum(reshape(vertex1(:,faces),[3 3 m]), 2));
    % Compute the Dirichlet energy of each face.
    E = zeros(1,m);
    for i=1:3
        i1 = mod(i,3)+1;
        % directed edge
        u = vertex1(:,faces(i,:)) - vertex1(:,faces(i1,:));
        % norm squared
        u = sum(u.^2);
        % weights between the vertices
        w = W(faces(i,:) + (faces(i1,:)-1)*n);
        E = E + w.*u;
    end
    % Compute gradient direction.
    G = zeros(3,n);
    Edir(end+1) = 0;
    for j=1:m
        f = faces(:,j);
        Alpha = A(:,j);
        alpha = norm(Alpha);
        for i=1:3
            i1 = mod(i  ,3)+1;
            i2 = mod(i+1,3)+1;
            % directed edges
            u1 = vertex1(:,f(i)) - vertex1(:,f(i1));
            u2 = vertex1(:,f(i)) - vertex1(:,f(i2));
            % weights between the vertices
            w1 = W(f(i) + (f(i1)-1)*n);
            w2 = W(f(i) + (f(i2)-1)*n);
            G(:,f(i)) = G(:,f(i)) + (w1*u1 + w2*u2) ./ alpha^2 - Alpha/alpha^4 * E(j);
        end
        Edir(end) = Edir(end) + E(j)/alpha;
    end
    % Perform the gradient descent step and the projection.    
    vertex1 = vertex1 - eta*G;
    vertex1 = vertex1 ./ repmat( sqrt(sum(vertex1.^2,1)), [3 1] );
    % display
    if 0
        if mod(it,5)==1
            clf;
            plot_mesh(vertex1,faces,options);
            axis tight; shading faceted;
            drawnow;
        end
        
    end
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
