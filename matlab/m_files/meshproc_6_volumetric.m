%% Volumetric Meshes
% This tour explores the processing of volumetric tetrahedral meshes.

perform_toolbox_installation('signal', 'general', 'graph', 'additional');

%% Tetrahedral Mesh Loading and Displaying
% You can load and display volumetric tetrahedral meshes.
% Important: .tet files and *not* included in the toolbox distribution (too
% large files). You should download them from 

%%
% http://www.aimatshape.net/

%%
% Load a volumetric mesh.

[vertex,faces] = read_tet('hand.tet');

%%
% Display it.

clear options;
options.plot_points = 1;
clf; plot_mesh(vertex,faces,options);

%%
% Display it.

options.cutting_plane = [0 0 1];
options.plot_points = 0;
clf; plot_mesh(vertex,faces,options);

%%
% Another view.

options.cutting_plane = [0 -1 0];
options.plot_points = 0;
options.cutting_offs = -.2;
options.face_vertex_color = vertex(1,:)';
clf; plot_mesh(vertex,faces,options);
view(-20,45);zoom(.8);
colormap jet(256);

%% Difusion Inside the Volume
% One can compute averaging operator inside the mesh.

%EXO
%% Compute the combinatorial adjacency matrix |W|.
% Compute edge list.
u = [1 2; 2 3; 3 1; 1 4; 2 4; 3 4]';
E = faces(u,:); E = reshape(E, [2 6*size(E,2)]);
p = size(E,2);
% Compute the adjacency matrix.
W = sparse( E(1,:), E(2,:), ones(p,1) );
W = max(W,W');
%EXO

%%
% Compute the combinatorial Laplacian, stored as a sparse matrix.

n = size(W,1);
D = spdiags(sum(W)', 0, n,n);
L = D-W;

%%
% Compute a set of random sources.

i = round( rand(20,1)*n ) + 1;
b = zeros(n,1); 
b(i) = (-1).^(1:length(i));

%% 
% Compute the diffusion with fixed boundary conditions.


L1 = L;
L1(i,:) = 0; L1(i+(i-1)*n) = 1;
v = L1\b;

%%
% Display.

options.face_vertex_color = v;
clf; plot_mesh(vertex,faces,options);
view(-20,45);zoom(.8);
shading interp;
colormap jet(256);
