niter = 12;
vertex1 = vertex;
clf;
for i=1:niter
    % Compute the delaunay triangulation.
    faces1 = delaunay(vertex1(1,:),vertex1(2,:))';
    % Compute the list of edges.
    E = [faces1([1 2],:) faces1([2 3],:) faces1([3 1],:)];
    p = size(E,2);
    % We build the adjacency matrix of the triangulation.
    A = sparse( E(1,:), E(2,:), ones(p,1) );
    % Normalize the adjacency matrix to obtain a smoothing operator.
    d = 1./sum(A);
    iD = spdiags(d(:), 0, n,n);
    W = iD * A;
    for q=1:1
    % Apply the filtering and ensure unit variance.
    vertex1 = vertex1*W';
    % force position
    vertex1(:,1:m) = vertexF;
    end
    % clf; plot_mesh(vertex1,faces1); drawnow;
    % Display the positions before / after.
    if mod(i,ceil(niter/4))==0
        subplot(2,2,i/ceil(niter/4));
        plot_mesh(vertex1,faces1);
        title(['Iteration ' num2str(i)]);
    end
end
