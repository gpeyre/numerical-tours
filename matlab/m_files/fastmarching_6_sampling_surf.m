%% Geodesic Surface Remeshing 
% This tour explores geodesic remeshing of surfaces.

%%
% This method is introduced in

%%
% _Geodesic Remeshing Using Front Propagation_
% Gabriel Peyré and Laurent Cohen, 
% International Journal on Computer Vision, Vol. 69(1), p.145-156, Aug. 2006.

perform_toolbox_installation('signal', 'general', 'graph');


%% Farthest Point Sampling
% An uniform sampling of points on a surface is obtained using a greedy
% farthest point sampling.

%%
% Load a 3D mesh.
    
clear options;
name = 'bunny';
[vertex,faces] = read_mesh(name);
n = size(vertex,2);
options.name = name;

%%
% Display it.

clf;
plot_mesh(vertex,faces, options);

%%
% Pick a first point.

landmarks = [100];

%%
% Compute the geodesic distance to this point.

[D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks);

%%
% Display the geodesic distance to the point.

clf; hold on;
options.face_vertex_color = mod( 20*D/max(D),1 );
plot_mesh(vertex,faces, options);
colormap jet(256);
h = plot3(vertex(1,landmarks), vertex(2,landmarks), vertex(3,landmarks), 'r.');
set(h, 'MarkerSize', 20);

%%
% Select as the next sampling point the farthest point.

[tmp,landmarks(end+1)] = max(D);

%%
% Update the distance map using a local propagation.

options.constraint_map = D;
[D1,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks,options);
D = min(D,D1);

%%
% Display the update distance map.

clf; hold on;
options.face_vertex_color = mod( 20*D/max(D),1 );
plot_mesh(vertex,faces, options);
colormap jet(256);
h = plot3(vertex(1,landmarks), vertex(2,landmarks), vertex(3,landmarks), 'r.');
set(h, 'MarkerSize', 20);

%EXO
%% Perform the farthest point sampling of |m=500| points.
m = 500;
% initialize
landmarks = [100];
[D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks);
clf;
k = 1; displist = [10 50 100 m];
for i=2:m
    % select
    [tmp,landmarks(end+1)] = max(D);
    % update
    options.constraint_map = D;
    [D1,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks,options);
    D = min(D,D1);
    if i==displist(k)
        subplot(2,2,k);
        hold on;
        options.face_vertex_color = perform_hist_eq(D,'linear');
        plot_mesh(vertex,faces, options);
        colormap jet(256);
        h = plot3(vertex(1,landmarks), vertex(2,landmarks), vertex(3,landmarks), 'r.');
        set(h, 'MarkerSize', 20);
        k = k+1;
    end
end
%EXO


%% Geodesic Delaunay Triangulation
% An intrinsic triangulation of the point is obtained using the geodesic
% Delaunay triangulation.

%%
% Compute the voronoi map |Q| of the segmentation.

[D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks);

%%
% Display the update distance map.

[B,I,J] = unique(Q);
v = randperm(m)'; J = v(J);
clf; hold on;
options.face_vertex_color = J;
plot_mesh(vertex,faces, options);
colormap jet(256);
h = plot3(vertex(1,landmarks), vertex(2,landmarks), vertex(3,landmarks), 'k.');
set(h, 'MarkerSize', 15);

%%
% Count the number |d(i)| of different voronoi indexes for each face |i|.

V = Q(faces); V = sort(V,1);
V = unique(V', 'rows')';
d = 1 + (V(1,:)~=V(2,:)) + (V(2,:)~=V(3,:));

%%
% Select the faces with 3 different indexe, they corresponds to geodesic
% Delaunay faces.

I = find(d==3); I = sort(I);

%%
% Build the Delaunay faces set.

z = zeros(n,1);
z(landmarks) = (1:m)';
facesV = z(V(:,I));

%%
% Position of the vertices of the subsampled mesh.

vertexV = vertex(:,landmarks);

%%
% Re-orient the faces so that they point outward of the mesh.

options.method = 'slow';
options.verb = 0;
facesV = perform_faces_reorientation(vertexV,facesV, options);

%%
% Display the sub-sampled mesh.

clf;
options.face_vertex_color = [];
plot_mesh(vertexV,facesV, options);
shading faceted;


%% Spacially Varying Remeshing
% It is possible to seed more point on a given part of the mesh.

%% 
% Create a density function by designing an isotropic metric.
% Here we use a metric that is slower in the left part.

W = ones(n,1);
W(vertex(1,:)<median(vertex(1,:))) = .4;
options.W = W;

%%
% Display the speed function.

clf;
hold on;
options.face_vertex_color = W;
plot_mesh(vertex,faces, options);
colormap jet(256);

%%
% Perform front propagation using this speed function.

landmarks = [5000];
options.constraint_map = [];
[D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks, options);

%%
% Display the distance map.

clf;
hold on;
options.face_vertex_color = mod( 20*D/max(D),1 );
plot_mesh(vertex,faces, options);
colormap jet(256);
h = plot3(vertex(1,landmarks), vertex(2,landmarks), vertex(3,landmarks), 'r.');
set(h, 'MarkerSize', 20);

%EXO
%% Perform a spacially adative remeshing.
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
%EXO


%% Feature Sensitive Remeshing
% A better remeshing quality is obtained by sampling more densly sharp
% features. This is achieved using a spatially varying metric, so that the front propagate slowly near
% regions of high curvature.

%%
% Compute the curvature of the mesh.

[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(vertex,faces,options);
 
%%
% Compute the total curvature.
 
C = abs(Cmin)+abs(Cmax);
 
%%
% Display it.

clf;
hold on;
options.face_vertex_color = min(C,.1);
plot_mesh(vertex,faces, options);
colormap jet(256);

%EXO 
%% Design a metric |W| so that the sampling is densed in area where |C| is
%% large.
W = rescale( min(C,.1), .001, 1);
options.W = W;
%
landmarks = [5000];
options.constraint_map = [];
[D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks, options);
% display
clf;
hold on;
options.face_vertex_color = mod( 20*D/max(D),1 );
plot_mesh(vertex,faces, options);
colormap jet(256);
h = plot3(vertex(1,landmarks), vertex(2,landmarks), vertex(3,landmarks), 'r.');
set(h, 'MarkerSize', 20);
%EXO

%EXO
%% Use such a metric to perform feature sensitive remeshing.
%% Tune the metric to reduce as much as possible the Hausdorff
%% approximation error.
%EXO