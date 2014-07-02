%% Spherical Haar Wavelets
% This tour explores multiscale computation on a 3D multiresolution sphere
% using a face-based haar transform.

perform_toolbox_installation('signal', 'general', 'graph', 'wavelet_meshes/');

%CMT
rep = 'results/mesh_haar_sphere/';
if not(exist(rep))
    mkdir(rep);
end
addpath('toolbox_additional/');
%CMT

%% Functions Defined on Surfaces
% One can define a function on a discrete 3D mesh that assigns a value to
% each vertex. One can then perform processing of the function according to
% the geometry of the surface. Here we use a simple sphere.


%%
% First compute a multiresolution sphere.

options.base_mesh = 'ico';
options.relaxation = 1;
options.keep_subdivision = 1;
J = 7;
[vertex,face] = compute_semiregular_sphere(J,options);
n = size(face{end},2);

%%
% Display two examples of sphere.

clf;
for j=[1 2 3 4]
    subplot(2, 2, j);
    plot_mesh(vertex{j}, face{j});
    shading faceted;
end

%%
% Comput the center of each face.

x = [];
for i=1:3
    v = vertex{end}(i,:);
    x(i,:) = mean(v(face{end}));
end


%%
% Load an image.

name = 'lena';
%CMT
name = 'haar';
%CMT
M = rescale( load_image(name, 512) );

%%
% Display it.

clf;
imageplot(crop(M));

%%
% Load a function on the sphere.
% Use the center of each face to sample the function.

f = rescale( load_spherical_function(M, x, options) );

%% 
% Display the function on the sphere.

vv = [125,-15];
%CMT
vv = [30 5];
%CMT
options.face_vertex_color = f;
clf;
plot_mesh(vertex{end}, face{end}, options);
view(vv);
colormap gray(256);
lighting none;

%CMT
mysaveas = @(a)saveas(gcf, [rep name '-' a '.png'], 'png');
myimwrite = @(a,x)imwrite(clamp(x), [rep name '-' a '.png'], 'png');
mysaveas('sphere-original');
%CMT

%% Multiscale Low Pass
% One can compute low pass approximation by iteratively averaging over 4
% neighboring triangles.

%%
% Perform one low pass filtering.

f1 = mean( reshape(f, [length(f)/4 4]), 2);

%%
% Display.

clf;
options.face_vertex_color = f1;
plot_mesh(vertex{end-1}, face{end-1}, options);
view(vv);
lighting none;

%EXO
%% Compute the successive low pass approximations.
fj = f;
clf;
for j=1:4
    fj = reshape(fj, [length(fj)/4 4]);
    fj = mean(fj,2);
    subplot(2,2,j);
    options.face_vertex_color = fj;
    plot_mesh(vertex{end-j}, face{end-j}, options);
    view(vv);
    colormap gray(256);
    lighting none;
end
%EXO

%CMT
fj = f;
for j=1:4
    fj = reshape(fj, [length(fj)/4 4]);
    fj = mean(fj,2);
    options.face_vertex_color = fj;
    clf;
    plot_mesh(vertex{end-j}, face{end-j}, options);
    view(vv);
    colormap gray(256);
    lighting none;
    mysaveas(['square-linapprox-' num2str(j)] );
end
%CMT


%% Spherical Haar Transform
% One can compute a wavelet transform by extracting, at each scale, 3
% orthogonal wavelet coefficient to represent the orthogonal complement
% between the successive resolutions.

%% 
% Precompute the local wavelet matrix, which contains the local vector and
% three orthognal detail directions.

U = randn(4);
U(:,1) = 1;
[U,R] = qr(U);

%%
% Initialize the forward transform.

fw = f;
nj = length(f);

%% 
% Extract the low pass component and apply the matrix U

fj = fw(1:nj);
fj = reshape(fj, [nj/4 4]);
fj = fj*U;

%%
% Store back the coefficients.

fw(1:nj) = fj(:);
nj = nj/4;    

%EXO 
%% Compute the full wavelet transform, and check for orthogonality
%% (conservation of energy).
fw = f;
nj = length(f);
for j=1:J-1
    fj = fw(1:nj);
    fj = reshape(fj, [nj/4 4]);
    fj = fj*U;
    fw(1:nj) = fj(:);
    nj = nj/4;
end
% Check for orthogonality.
e = abs(norm(f)-norm(fw))/norm(f);
disp(['Orthogonality deviation (should be 0): ' num2str(e)]);
%EXO

%%
% Display the coefficients "in place".

clf;
options.face_vertex_color = clamp(fw,-2,2);
plot_mesh(vertex{end}, face{end}, options);
view(vv);
colormap gray(256);
lighting none;

%CMT
mysaveas('sphere-wavelet');
%CMT

%%
% Display the decay of the coefficients.

clf;
plot(fw);
axis([1 n -5 5]);

%EXO
%% Implement the backward spherical Haar transform (replace U by U' to perform the reconstruction), 
%% and check for perfect reconstruction.
nj = size(face{1},2);
f1 = fw;
for j=1:J-1
    fj = f1(1:4*nj);
    fj = reshape(fj, [nj 4]);
    fj = fj*U';
    f1(1:4*nj) = fj(:);
    nj = nj*4;
end
e = norm(f-f1)/norm(f);
disp(['Bijectivity deviation (should be 0): ' num2str(e)]);
%EXO


clf;
options.face_vertex_color = clamp(f1);
plot_mesh(vertex{end}, face{end}, options);
view(vv);
colormap gray(256);
lighting none;

%EXO
%% Perform Haar wavelet approximation with only 10% of the coefficients.
% Forward transform
rho = .1;
rho = .05;
nj = length(f);
fw = f;
for j=1:J-1
    fj = fw(1:nj);
    fj = reshape(fj, [nj/4 4]);
    fj = fj*U;
    fw(1:nj) = fj(:);
    nj = nj/4;
end
% thresholding
fwT = perform_thresholding(fw, round(rho*n), 'strict');
% Backward transform
nj = size(face{1},2);
f1 = fwT;
for j=1:J-1
    fj = f1(1:4*nj);
    fj = reshape(fj, [nj 4]);
    fj = fj*U';
    f1(1:4*nj) = fj(:);
    nj = nj*4;
end
% display
clf;
options.face_vertex_color = clamp(f1);
plot_mesh(vertex{end}, face{end}, options);
view(vv);
colormap gray(256);
lighting none;
%EXO

%CMT
mysaveas(['sphere-approx-' num2str(round(100*rho))]);
%CMT

%EXO
%% Compare with the traditional 2D Haar approximation of |M|.
options.h = compute_wavelet_filter('Daubechies',2);
MW = perform_wavortho_transf(M,0,+1,options);
MWT = perform_thresholding(MW(:), round(rho*n), 'strict');
MWT = reshape(MWT, size(M));
M1 = perform_wavortho_transf(MWT,0,-1,options);
clf;
imageplot(clamp(crop(M1)));
%EXO

%CMT
myimwrite('square-original', crop(M, 256, [.5 .3]*512) );
myimwrite(['square-approx-' num2str(round(100*rho))], crop(M1, 256, [.5 .3]*512) );
%CMT

%EXO
%% Implement Spherical denoising using the Haar transform. Compare it with
%% vertex-based lifting scheme denoising.
%EXO