%% Geodesic Triangulation for Image Compression
% This tour explores the use geodesic triangulations to perform image
% compression.

perform_toolbox_installation('signal', 'general', 'graph');

%% Image Approximation with Triangulation
% It is possible to approximate an image over a triangulation using
% piecewise linear splines. 

%%
% Load an image.

name = 'cameraman';
n = 256;
M = rescale( load_image(name, n) );

%%
% Number of points used to compute the approximation. The more points, the
% smallest the error.

m = 400;

%%
% Seeds random points, include the corners into these points.

vertex = floor( rand(2,m-4)*(n-1) ) +1;
vertex(:,end+1:end+4) = [[1;1] [1;n] [n;n] [n;1]];

%%
% Compute a Delaunay triangulation.

faces = compute_delaunay(vertex);

%%
% A first way to perform an approximation with |m| triangles is to
% interpolate the image at the sampling points |vertex|.

vinterp = interp2(M, vertex(2,:), vertex(1,:));

%%
% Each |vinterp(i)| is the value of the approximating image at the 
% |vertex(:,i)|. We compute a spline interpolation.

Minterp = compute_triangulation_interpolation(faces,vertex,vinterp, n);

%%
% Display the approximation.

clf;
subplot(1,2,1);
plot_triangulation(vertex,faces, M);
title('Triangulation');
subplot(1,2,2);
imageplot(clamp(Minterp), ['Interpolation, SNR=' num2str(snr(Minterp,M)) 'dB']);

%%
% Another, better way to compute the approximation is to compute
% coefficients |vapprox| that performs the best L2 approximation with
% linear spline.

vapprox = compute_orthoproj_triangulation(vertex, faces, M);
Mapprox = compute_triangulation_interpolation(faces,vertex,vapprox, n);

%%
% Compare interpolation and approximation.

clf;
imageplot(clamp(Minterp), ['Interpolation, SNR=' num2str(snr(Minterp,M),3) 'dB'], 1,2,1);
imageplot(clamp(Mapprox), ['Approximation, SNR=' num2str(snr(Mapprox,M),3) 'dB'], 1,2,2);

%% Isotropic Metrics for Image Approximation
% It is possible to compute optimized sampling location |vertex| by using
% the farthest point sampling algorithm with a well chosen metric |W| so
% that more points are put in areas of strong gradient.

%%
% The metric will be of the form

%%
% |W(x) = (norm(grad_x(M))+epsilon|)^alpha|

%%
% where |epsilon| and |alpha| control the density variation strength.

%%
% Parameters for the metric.

alpha = .7;
epsilon = 1e-2;

%EXO
%% Compute a density function that is larger at area of large gradient.
%% |W(x) = (norm(grad(M))+epsilon|)^alpha|, for |alpha=.7|.
%% To stabilize the process, you can smooth a bit the gradient magnitude.
options.order = 2;
G = grad(M, options);
W = sum(G.^2,3);
% blur a little
W = perform_blurring(W, 10);
W = (W+epsilon).^alpha;
% scale to set up the contast
clf;
imageplot(M, 'Image', 1,2,1);
imageplot(W, 'Metric', 1,2,2);
%EXO

%EXO
%% Perform farthest points sampling to compute sampling location |vertex| 
%% and the corresponding geodesic Delaunay triangulation |faces|.
ndisp = [20 50 100 m];
vertex = [1;1];
[D,Z,Q] = perform_fast_marching(1./W, vertex);
vertex_svg = {};
faces_svg = {};
clf; u = 1;
for k=2:m
    options.constraint_map = D;
    [D1,Z,Q] = perform_fast_marching(1./W, vertex(:,end), options);
    D = min(D,D1);
    [tmp,i] = max(D(:));
    [x,y] = ind2sub([n n],i); 
    vertex(:,end+1) = [x;y];
    if k==ndisp(u)
        subplot(2,2,u);
        hold on;
        imageplot(M, [num2str(k) ' points']);
        plot(vertex(2,:), vertex(1,:), 'r.');
        axis ij;
        u = u+1;  
    end
end
% compute the Delaunay triangulation
[D,Z,Q] = perform_fast_marching(1./W, vertex);
faces = compute_voronoi_triangulation(Q,vertex);
%EXO

%%
% Perform approximation over the triangulation.

vgeod = compute_orthoproj_triangulation(vertex, faces, M);
Mgeod = compute_triangulation_interpolation(faces,vertex,vgeod, n);

%%
% Compare interpolation and approximation.

clf;
subplot(1,2,1);
plot_triangulation(vertex,faces, M);
subplot(1,2,2);
imageplot(clamp(Mgeod), ['SNR=' num2str(snr(Mgeod,M),3) 'dB']);

%EXO
%% For a large value of |m| compute the approximation for several |alpha|.
alpha_list = [.4 .7 1 1.3];
m = 600;
clf;
for ialpha=1:length(alpha_list)
    alpha = alpha_list(ialpha);
    % metric
    W = sum(G.^2,3);
    W = perform_blurring(W, 6);
    W = (W+epsilon).^alpha;
    % Farthest point
    vertex = [1;1];
    [D,Z,Q] = perform_fast_marching(1./W, vertex);
    for k=2:m
        options.constraint_map = D;
        [D1,Z,Q] = perform_fast_marching(1./W, vertex(:,end), options);
        D = min(D,D1);
        [tmp,i] = max(D(:));
        [x,y] = ind2sub([n n],i);
        vertex(:,end+1) = [x;y];
    end
    % compute the Delaunay triangulation
    [D,Z,Q] = perform_fast_marching(1./W, vertex);
    faces = compute_voronoi_triangulation(Q,vertex);
    % compute approximation
    vgeod = compute_orthoproj_triangulation(vertex, faces, M);
    Mgeod = compute_triangulation_interpolation(faces,vertex,vgeod, n);
    % display
    imageplot(clamp(Mgeod), ['\alpha= ' num2str(alpha) ', SNR=' num2str(snr(Mgeod,M),3) 'dB'], 2,2,ialpha);
end
%EXO
