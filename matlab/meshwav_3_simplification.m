%% Mesh Simplification
% This tour explore the simplication of a highly detailed mesh into a
% coarser one.

perform_toolbox_installation('signal', 'general', 'graph');

%% Random Edge Contraction
% Simplest way to perform mesh simplification is through edge collapse.
% Each step replaces two vertex joined by an edge by a single vertex, and
% removes the edge.

%%
% Load a model.

name = 'venus';
options.name = name;
[vertex,faces] = read_mesh(name);
n = size(vertex,2);

%%
% Display full quality.

plot_mesh(vertex,faces,options);
shading faceted;

%%
% Initialize the simplification.

faces1 = faces;
vertex1 = vertex;

%% 
% Compute the collection of edges.

edges = compute_edges(faces1);
nedges = size(edges,2);

%% 
% Select an edge. Here selection is done at random.

k = floor(rand*(nedges-1))+1;
e = edges(:,k);

%%
% Change the vertex location, and remove one of the two vertices.

vertex1(:,e(1)) = mean( vertex1(:,e),2 );
vertex1(:,e(2)) = Inf;

%%
% Change the face indexing.

faces1(faces1==e(2)) = e(1);
a = sum( diff(sort(faces1))==0 );
faces1(:,a>0) = [];

    
%EXO
%% Perform iterative collapse to reach |p = round(2*n/3)| vertices.
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
%EXO

%EXO
%% As a post processing, find a way to remove from |faces1| and |vertex1| the unecessary
%% information (remove vertex and faces that are not used).
%EXO


%% Error Driven Edge Contraction
% To ameliorate the quality of the output mesh, it is necessary to order to
% select the edge to collapse according to some quality priority.

%EXO
%% Perform iterative collapse to reach |p = round(2*n/3)| vertices.
%% Use an ordering of the edge by their length.
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
%EXO

%EXO 
%% Try to use other criteria.
% No correction for this exercise.
%EXO