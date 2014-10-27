%% Geodesic Medial Axsis
% This tour studies the computation of the medial axis using the Fast
% Marching.

perform_toolbox_installation('signal', 'general', 'graph');


%% Voronoi Diagram
% The Voronoi diagram is the segmentation of the image given by the region
% of influence of the set of starting points.

%%
% Load a distance map.

n = 200;
W = load_image('mountain', n);
W = rescale(W,.25,1);

%% 
% Select seed points.

pstart = [[20;20] [120;100] [180;30] [60;160]];
nbound = size(pstart,2);

%%
% Display the map and the points.

ms = 20;
clf; hold on;
imageplot(W);
h = plot(pstart(2,:), pstart(1,:), '.r'); set(h, 'MarkerSize', ms);

%%
% Compute the geodesic distant to the whole set of points.

[D,S,Q] = perform_fast_marching(W, pstart);

%%
% Display the geodesic distance.

clf; hold on;
imageplot(convert_distance_color(D, W));
h = plot(pstart(2,:), pstart(1,:), '.r'); set(h, 'MarkerSize', ms);

%%
% Display the Voronoi Segmentation.

clf; hold on;
imageplot(Q);
h = plot(pstart(2,:), pstart(1,:), '.r'); set(h, 'MarkerSize', ms);
colormap jet(256);

%% Medial Axis from the Voronoi Map
% The medial axis is difficult to extract from the singularity of the
% distance map. It is much more robust to extract it from the
% discontinuities in the Voronoi index map |Q|.

%%
% Compute the derivative, the gradient.

G = grad(Q);

%%
% Take it modulo |nbound|.

G(G<-nbound/2) = G(G<-nbound/2) + nbound;
G(G>nbound/2) = G(G>nbound/2) - nbound;

%%
% Compute the norm of the gadient.

G = sqrt(sum(G.^2,3));

%%
% Compute the medial axis by thresholding the gradient magnitude.

B = 1 - (G>.1);

%%
% Display.

clf; hold on;
imageplot(B);
h = plot(pstart(2,:), pstart(1,:), '.r'); set(h, 'MarkerSize', ms);



%% Skeleton of a Shape 
% The sekeleton, also called Medial Axis, is the set of points where the
% geodesic distance is singular. 

%%
% A binary shape is represented as a binary image.

n = 200;
name = 'chicken';
M = load_image(name,n);
M = perform_blurring(M,5);
M = double( rescale( M )>.5 );
if M(1)==1 
    M = 1-M;
end

%%
% Compute its boundary, that is going to be the set of starting points.

pstart = compute_shape_boundary(M);
nbound = size(pstart,2);

%%
% Display the metric.

lw = 2;
clf; hold on;
imageplot(-M);
h = plot(pstart(2,:), pstart(1,:), 'r'); set(h, 'LineWidth', lw); axis ij;

%% 
% Parameters for the Fast Marching: constant speed |W|, but retricted using |L| to the
% inside of the shape.

W = ones(n);
L = zeros(n)-Inf; L(M==1) = +Inf;


%%
% Compute the fast marching, from the boundary points.

options.constraint_map = L;
[D,S,Q] = perform_fast_marching(W, pstart, options);
D(M==0) = Inf;

%%
% Display the distance function to the boundary.

clf;
hold on;
display_shape_function(D);
h = plot(pstart(2,:), pstart(1,:), 'r'); set(h, 'LineWidth', lw); axis ij;


%%
% Display the index of the closest boundary point.

clf;
hold on;
display_shape_function(Q);
h = plot(pstart(2,:), pstart(1,:), 'r'); set(h, 'LineWidth', lw); axis ij;

%EXO
%% Compute the norm of the gradient |G| modulo |nbound|. Be careful to remove the boundary
%% of the shape from this indicator. Display the thresholded gradient map.
% gradient
G = grad(Q);
G(G<-nbound/2) = G(G<-nbound/2) + nbound;
G(G>nbound/2) = G(G>nbound/2) - nbound;
% Compute the norm of the gadient.
G = sqrt(sum(G.^2,3));
% Remove the boundary to the skeletton.
M1 = perform_convolution(M,ones(3)/9)>.99;
G = G.*M1;
clf;
B = G>15;
A = M;
A(B==1) = 0;
imageplot(-A);
%EXO


%EXO
%% Display the Skeleton obtained for different threshold values.
thresh = [10 20 50 100];
clf;
for i=1:4
    subplot(2,2,i);
    B = G>thresh(i);
    A = M;
    A(B==1) = 0;
    imageplot(-A);
end
%EXO

%% Skeletton