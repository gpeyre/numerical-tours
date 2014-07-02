%% Basics About 2D Triangulation
% This tour explores some basics about 2D triangulated mesh (loading,
% display, manipulations).
perform_toolbox_installation('signal', 'general', 'graph');


%% Planar Triangulation
% A planar triangulation is a collection of |n| 2D points, whose coordinates are
% stored in a |(2,n)| matrix |vertex|, and a topological collection of triangle, stored
% in a |(m,2)| matrix |faces|.

%%
% Number of points.

n = 200;

%%
% Compute randomized points in a square.

vertex = 2*rand(2,n)-1;

%%
% A simple way to build a triangulation of the convex hull of the points is
% to compute the Delaunay triangulation of the points.

faces = delaunay(vertex(1,:),vertex(2,:))';

%%
% One can display the triangulation.

clf;
subplot(1,2,1);
hh = plot(vertex(1,:),vertex(2,:), 'k.');
axis('equal'); axis('off');
set(hh,'MarkerSize',10);
title('Points');
subplot(1,2,2);
plot_mesh(vertex,faces);
title('Triangulation');

%% Point Modification
% It is possible to modify the position of the points like a particles
% system. The dynamics is govered by the connectivity to enfoce an even distribution.
% During the modification of the positions, the connectivity is updated.

%%
% Fix some points on a disk.

m = 20;
t = linspace(0,2*pi,m+1); t(end) = [];
vertexF = [cos(t);sin(t)];
vertex(:,1:m) = vertexF;
faces = delaunay(vertex(1,:),vertex(2,:))';

%%
% Initialize the positions.

vertex1 = vertex;

%%
% Compute the delaunay triangulation.

faces1 = delaunay(vertex1(1,:),vertex1(2,:))';

%% 
% Compute the list of edges.

E = [faces([1 2],:) faces([2 3],:) faces([3 1],:)];
p = size(E,2);

%%
% We build the adjacency matrix of the triangulation.

A = sparse( E(1,:), E(2,:), ones(p,1) );

%%
% Normalize the adjacency matrix to obtain a smoothing operator.

d = 1./sum(A);
iD = spdiags(d(:), 0, n,n);
W = iD * A;

%%
% Apply the filtering.

vertex1 = vertex1*W';

%%
% Set of the position of fixed points.

vertex1(:,1:m) = vertexF;


%%
% Display the positions before / after.

clf;
subplot(1,2,1);
plot_mesh(vertex,faces);
title('Before filering');
subplot(1,2,2);
plot_mesh(vertex1,faces1);
title('After filtering');

%EXO
%% Compute several steps of iterative filterings, while ensuring the positions of the fixed points.
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
%EXO
