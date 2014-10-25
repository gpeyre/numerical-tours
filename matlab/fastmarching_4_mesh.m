%% Geodesic Mesh Processing 
% This tour explores geodesic computations on 3D meshes.

perform_toolbox_installation('signal', 'general', 'graph', 'wavelet_meshes');


%% Distance Computation on 3D Meshes
% Using the fast marching on a triangulated surface, one can compute the
% distance from a set of input points. 
% This function also returns the segmentation of the surface into geodesic Voronoi cells.

%%
% Load a 3D mesh.
    
name = 'elephant-50kv';
[vertex,faces] = read_mesh(name);
nvert = size(vertex,2);

%%
% Starting points for the distance computation.

nstart = 15;
pstarts = floor(rand(nstart,1)*nvert)+1;
options.start_points = pstarts;

%%
% No end point for the propagation.

clear options;
options.end_points = [];

%%
% Use a uniform, constant, metric for the propagation.

options.W = ones(nvert,1);

%%
% Compute the distance using Fast Marching.

options.nb_iter_max = Inf;
[D,S,Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options);

%% 
% Display the distance on the 3D mesh.

clf;
plot_fast_marching_mesh(vertex,faces, D, [], options);

%% 
% Extract precisely the voronoi regions, and display it.

[Qexact,DQ, voronoi_edges] = compute_voronoi_mesh(vertex, faces, pstarts, options);
options.voronoi_edges = voronoi_edges;
plot_fast_marching_mesh(vertex,faces, D, [], options); 

%EXO
%% Using |options.nb_iter_max|, display the progression of the propagation.
clf;
nblist = round( linspace(.1,1,6)*nvert );
for i = 1:length(nblist);
    options.nb_iter_max = nblist(i);
    [D,S,Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options);
    subplot(2,3,i);
    col = D; col(col==Inf) = 0;
    options.face_vertex_color = col;
    hold('on');
    plot_mesh(vertex,faces,options);
    colormap jet(256);
    % display here starting points
end
%EXO

%% Geodesic Path Extraction
% Geodesic path are extracted using gradient descent of the distance map.

%%
% Select random endding points, from which the geodesic curves start.

nend = 40;
pend = floor(rand(nend,1)*nvert)+1;

%%
% Compute the vertices 1-ring.

vring = compute_vertex_ring(faces);

%EXO
%% For each point |pend(k)|, compute a discrete geodesic path |path| such
%% that |path(1)=pend(k)| and |D(path(i+1))<D(path(i))|
%% with |[path(i), path(i+1)]| being an edge of the mesh.
%% This means that |path(i+1)| is an element of |vring{path(i)}|.
%% Display the paths on the mesh.
pathsD = {};
for k=1:nend
    % path purely on edges
    vlist = pend(k);
    vprev = D(vlist(end));
    while true
        x0 = vlist(end);
        r = vring{x0};
        [v,J] = min(D(r));
        x = r(J);
        if v>=vprev || v==0
            break;
        end
        vprev = v;
        vlist(end+1) = x;
    end
    pathsD{end+1} = vertex(:,vlist);
end
plot_fast_marching_mesh(vertex,faces, Q, pathsD, options); 
%EXO

%%
% In order to extract a smooth path, one needs to use a gradient descent.

options.method = 'continuous';
paths = compute_geodesic_mesh(D, vertex, faces, pend, options);

%%
% Display the smooth paths.

plot_fast_marching_mesh(vertex,faces, Q, paths, options); 

%% Curvature-based Speed Functions
% In order to extract salient features of a surface, one can define a speed
% function that depends on some curvature measure of the surface.

%%
% Load a mesh with sharp features.

clear options;
name = 'fandisk';
[vertex,faces] = read_mesh(name);
options.name = name;
nvert = size(vertex,2);

%%
% Display it.

clf;
plot_mesh(vertex,faces, options);

%%
% Compute the curvature.

options.verb = 0;
[Umin,Umax,Cmin,Cmax] = compute_curvature(vertex,faces,options);

%%
% Compute some absolute measure of curvature.

C = abs(Cmin)+abs(Cmax);
C = min(C,.1);

%% 
% Display the curvature on the surface

options.face_vertex_color = rescale(C);
clf;
plot_mesh(vertex,faces,options);
colormap jet(256);

%%
% Compute a metric that depends on the curvature.
% Should be small in area that the geodesic should follow.

epsilon = .5;
W = rescale(-min(C,0.1), .1,1);

%% 
% Display the metric on the surface

options.face_vertex_color = rescale(W);
clf;
plot_mesh(vertex,faces,options);
colormap jet(256);

%%
% Starting points.

pstarts = [2564; 16103; 15840];
options.start_points = pstarts;


%%
% Compute the distance using Fast Marching.

options.W = W;
options.nb_iter_max = Inf;
[D,S,Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options);

%% 
% Display the distance on the 3D mesh.

options.colorfx = 'equalize';
clf;
plot_fast_marching_mesh(vertex,faces, D, [], options);

%EXO
%% Using |options.nb_iter_max|, display the progression of the propagation for constant |W|.
clf;
options.W = ones(nvert,1);
nblist = round( [.05 .15 .4 .6]*nvert );
for i = 1:length(nblist);
    options.nb_iter_max = nblist(i);
    [D0,S,Q0] = perform_fast_marching_mesh(vertex, faces, pstarts, options);
    subplot(2,2,i);
    plot_fast_marching_mesh(vertex,faces, D0, [], options);
end
%EXO

%EXO
%% Using |options.nb_iter_max|, display the progression of the propagation for a curvature based |W|.
clf;
options.W = W;
for i = 1:length(nblist);
    options.nb_iter_max = nblist(i);
    [D,S,Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options);
    subplot(2,2,i);
    plot_fast_marching_mesh(vertex,faces, D, [], options);
end
%EXO

%EXO
%% Extract geodesics.
% compute distances
options.nb_iter_max = Inf;
options.W = ones(nvert,1);
[D0,S,Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options);
options.W = W;
[D,S,Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options);
% compute paths
pend = [7175; 17991; 2293];
options.method = 'continuous';
paths0 = compute_geodesic_mesh(D0, vertex, faces, pend, options);
paths = compute_geodesic_mesh(D, vertex, faces, pend, options);
% display
clf;
subplot(1,2,1);
plot_fast_marching_mesh(vertex,faces, perform_hist_eq(D0, 'linear'), paths0, options); 
subplot(1,2,2);
plot_fast_marching_mesh(vertex,faces, perform_hist_eq(D, 'linear'), paths, options); 
%EXO

%% Texture-based Speed Functions
% One can take into account a texture to design the speed function.

clear options;
options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 0;
[vertex,faces] = compute_semiregular_sphere(7,options);
nvert = size(vertex,2);

%% 
% Load a function on the mesh.

name = 'earth';
f = load_spherical_function(name, vertex, options);
options.name = name;


%%
% Starting points.

pstarts = [2844; 5777];
options.start_points = pstarts;

%%
% Display the function.

clf;
plot_fast_marching_mesh(vertex,faces, f, [], options);
colormap gray(256);

%%
% Load and display the gradient magnitude of the function.

g = load_spherical_function('earth-grad', vertex, options);

%%
% Display it.

clf;
plot_fast_marching_mesh(vertex,faces, g, [], options);
colormap gray(256);

%%
% Design a metric.

W = rescale(-min(g,10),0.01,1);

%%
% Display it.

clf;
plot_fast_marching_mesh(vertex,faces, W, [], options);
colormap gray(256);


%EXO
%% Using |options.nb_iter_max|, display the progression of the propagation for a curvature based |W|.
clf;
options.W = W;
nblist = round( [.05 .15 .4 .6]*nvert );
for i = 1:length(nblist);
    options.nb_iter_max = nblist(i);
    [D,S,Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options);
    subplot(2,2,i);
    plot_fast_marching_mesh(vertex,faces, D, [], options);
end
%EXO

%EXO
%% Extract geodesics.
% compute distances
options.nb_iter_max = Inf;
options.W = ones(nvert,1);
[D0,S,Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options);
options.W = W;
[D,S,Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options);
% compute paths
pend = [9969; 5073];
options.method = 'continuous';
paths0 = compute_geodesic_mesh(D0, vertex, faces, pend, options);
paths = compute_geodesic_mesh(D, vertex, faces, pend, options);
% display
clf;
subplot(1,2,1);
plot_fast_marching_mesh(vertex,faces, f, paths0, options); 
subplot(1,2,2);
plot_fast_marching_mesh(vertex,faces, f, paths, options); 
colormap gray(256);
%EXO