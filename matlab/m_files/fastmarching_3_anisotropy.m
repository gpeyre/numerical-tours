%% Anisotropic Fast Marching
% This tour explores the use of geodesic distances for anisotropic metric.

perform_toolbox_installation('signal', 'general', 'graph');

%% Structure Tensor Field
% An anisotropic metric is given through a tensfor field |T| which
% is an |(n,n,2,2)| array, where |T(i,j,:,:)| is a positive definite
% symmetric matrix that define the metric at pixel |(i,j)|.

%%
% Here we extract the tensor field whose main eigenvector field is alligned
% with the direction of the texture. This can be achieved using the
% structure tensor field, which remove the sign ambiguity by tensorizing
% the gradient, and remove noise by filtering.

%%
% Load an image that contains an oscillating texture.

name = 'fingerprint';
n = 150;
M = rescale(load_image(name,n));

%%
% Compute its gradient.

options.order = 2;
G = grad(M,options);

%%
% Compute a rank-1 tensor with main eigenvector aligned in the direction
% orthogonal to the gradient.

T = zeros(n,n,2,2);
T(:,:,1,1) = G(:,:,2).^2;
T(:,:,2,2) = G(:,:,1).^2;
T(:,:,1,2) = -G(:,:,1).*G(:,:,2);
T(:,:,2,1) = -G(:,:,1).*G(:,:,2);

%%
% Smooth the field (blur each entry).

sigma = 12;
T = perform_blurring(T,sigma);

%%
% Compute the eigenvector and eigenvalues of the tensor field.

[e1,e2,l1,l2] = perform_tensor_decomp(T);

%%
% Display the main orientation field.

clf;
plot_vf(e1(1:10:n,1:10:n,:),M);
colormap(gray(256));

%% Anisotropic Fast Marching
% One can compute geodesic distance and geodesics using an anisotropic Fast
% Marching propagation.

%%
% Anisotropy of the tensor field.

anisotropy = .1;

%%
% Build the Riemannian metric using the structure tensor direction.

H = perform_tensor_recomp(e1,e2, ones(n),ones(n)*1/anisotropy );

%%
% Starting point.

pstart = [n n]/4;

%%
% Perform the propagation.

hx = 1/n; hy = 1/n;
H(i,j,:,:) = [1,0;0,1] * 10000000;
[D, dUx, dUy, Vor, L] = fm2dAniso([hx;hy], H, pstart);
D(1,1) = 0;

%%
% Display the result.

clf;
subplot(1,2,1);
imageplot(M, 'Image');
subplot(1,2,2);
hold on;
imageplot(convert_distance_color(D), 'Geodesic distance');
hh = plot(pstart(2),pstart(1), 'r.');
set(hh, 'MarkerSize',15);
axis('ij');
colormap(gray(256));

%EXO
%% Compute the geodesic distance for several anisotropy, and for several
%% starting points.
npoints = 10;
nstart = floor( rand(2,npoints)*n ) + 1;
aniso_list = [.02 .1 .2 .5];
clf;
for i=1:length(aniso_list)
    anisotropy = aniso_list(i);
    H = perform_tensor_recomp(e1,e2, ones(n),ones(n)*1/anisotropy );
    [D, dUx, dUy, Vor, L] = fm2dAniso([hx;hy], H, pstart);
    imageplot(convert_distance_color(D), ['Anisotropy=' num2str(anisotropy)],2,2,i);
end
%EXO

%% Farthest Point Sampling
% It is possible to perform an anisotropic sampling of the image using a
% Farthest point sampling strategy.

%%
% We use a highly anisotropic metric.

anisotropy = .02;
H = perform_tensor_recomp(e1,e2, ones(n),ones(n)*1/anisotropy );
    

%%
% Choose the initial point.

vertex = [1;1];

%%
% Compute the geodesic distance.

[D, dUx, dUy, Vor, L] = fm2dAniso([hx;hy], H, vertex);

%%
% Choose the second point as the farthest point.

[tmp,i] = max(D(:));
[x,y] = ind2sub([n n],i); 
vertex(:,end+1) = [x;y];

%%
% Display distance and points.

clf;
subplot(1,2,1);
hold on;
imageplot(M, 'Image'); axis ij;
plot(vertex(2,1), vertex(1,1), 'r.');
plot(vertex(2,2), vertex(1,2), 'b.');
subplot(1,2,2);
hold on;
imageplot( convert_distance_color(D), 'Distance'); axis ij;
plot(vertex(2,1), vertex(1,1), 'r.');
plot(vertex(2,2), vertex(1,2), 'b.');
colormap gray(256);


%%
% Update the value of the distance map with a partial propagation from the last added point.

[D1, dUx, dUy, Vor, L] = fm2dAniso([hx;hy], H, vertex);

%% 
% Display old/new.

clf;
imageplot( D, 'Old distance', 1,2,1 );
imageplot( D1, 'New distance', 1,2,2 );
colormap jet(256);

%% 
% Update.

D = D1;

%EXO
%% Perform farthest point sampling.
nmax = 200;
ndisp = [20 50 100 200];
vertex = [1;1];
[D, dUx, dUy, Vor, L] = fm2dAniso([hx;hy], H, vertex);
clf; u = 1;
vertex_svg = {};
faces_svg = {};
for k=2:nmax
    options.constraint_map = D;
    [D, dUx, dUy, Vor, L] = fm2dAniso([hx;hy], H, vertex);
    [tmp,i] = max(D(:));
    [x,y] = ind2sub([n n],i); 
    vertex(:,end+1) = [x;y];
    if k==ndisp(u)
        subplot(2,2,u);
        hold on;
        imageplot(D, [num2str(k) ' points']);
        plot(vertex(2,:), vertex(1,:), 'r.');
        u = u+1;  
        % compute the Delaunay triangulation
        % [D1,Z,Q] = perform_fast_marching(1./W, vertex);
        % vertex_svg{end+1} = vertex;
        % faces_svg{end+1} = compute_voronoi_triangulation(Q,vertex);
    end
end
%EXO

%% Anisotropic Image Approximation
% One can combine the farthest point sampling strategy with a geodesic
% Delaunay triangulation to perform image approximation.

%%
% Load an image.

n = 256;
name = 'peppers-bw';
M = rescale(load_image(name, n));

%EXO
%% Compute a metric |H| adapted to the approximation of this image.
%EXO

%EXO
%% Perform farthest point sampling.
%EXO

%EXO
%% Compute the geodesic Delaunay triangulation of this set of point.
%EXO

%EXO
%% Perform image approximation using linear splines.
%EXO