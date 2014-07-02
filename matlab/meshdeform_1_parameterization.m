%% Mesh Parameterization 
% This tour explores 2-D parameterization of 3D surfaces using linear methods.

%%
% A  review paper for mesh parameterization can be found in:

%%
% M.S. Floater and K. Hormann, 
% _Surface Parameterization: a Tutorial and Survey_
% in Advances in multiresolution for geometric modelling, p. 157-186, 2005.

%%
% See also:

%%
% K. Hormann,  K. Polthier and A. Sheffer
% _Mesh parameterization: theory and practice,_
% Siggraph Asia Course Notes

perform_toolbox_installation('signal', 'general', 'graph');

%CMT
rep = '../results/mesh_parameterization/';
if not(isdir(rep))
    mkdir(rep);
end
%CMT

%% Conformal Laplacian
% The conformal Laplacian uses the cotan weights to obtain an accurate
% discretization of the Laplace Beltrami Laplacian.

%%
% They where first introduces as a linear finite element approximation of the
% Laplace-Beltrami operator in:

%%
% U. Pinkall and K. Polthier, 
% _Computing discrete minimal surfaces and their conjugates_
% Experimental Mathematics, 2(1):15-36, 1993.

%%
% First load a mesh. The faces are stored in a matrix \(F = (f_j)_{j=1}^m \in \RR^{3 \times m}\)
% of \(m\) faces \(f_j \in \{1,\ldots,n\}^3\). The position of the vertices are
% stored in a matrix \(X  = (x_i)_{i=1}^n \in \RR^{3 \times n}\) of \(n\) triplets of points
% \(x_k \in \RR^3\)

name = 'nefertiti';
[X,F] = read_mesh(name);
n = size(X,2);
clear options; options.name = name;

%%
% In order to perform mesh parameterization, it is important that this mesh
% has the topology of a disk, i.e. it should have a single B.

%%
% First we compute the boundary \(B = (i_1,\ldots,i_p)\) of the mesh.
% By definition, for the edges \((i_k,i_{k+1})\), there is a single adjacent face \(
% (i_k,i_{k+1},\ell) \).

options.verb = 0;
B = compute_boundary(F, options);

%%
% Length of the boundary.

p = length(B);

%%
% Display the boundary.

clf; hold on;
plot_mesh(X,F,options);
shading('faceted');
hh = plot3( X(1,B),X(2,B),X(3,B), 'r' );
set_linewidth(hh,4);


%%
% The conformal Laplacian weight matrix \(W \in \RR^{n \times n}\) is
% defined as 
% \[ W_{i,j} = 
%       \choice{
%           \text{cotan}(\al_{i,j}) + \text{cotan}(\be_{i,j})
%           \qifq i \sim j \\
%           \quad 0 \quad \text{otherwise}.
%       }
% \]
% Here, \(i \times j\) means that there exists two faces
% \( (i,j,k) \) and \( (i,j,\ell) \) in the mesh (note that for B
% faces, one has \(k=\ell\)).

%%
% The angles are the angles centered as \(k\) and \(\ell\), i.e.
% \[ \al_{i,j} = \widehat{x_i x_k x_j }
%   \qandq
%   \be_{i,j} = \widehat{x_i x_\ell x_j }. \]

%%
% Compute the conformal 'cotan' weights. 
% Note that each angle \(\alpha\) in the mesh contributes with 
% \(1/\text{tan}(\alpha)\) to the weight of the opposite edge.
% We compute \(\alpha\) as
% \[ \alpha = \text{acos}\pa{ \frac{\dotp{u}{v}}{\norm{u}\norm{v}} } \]
% where \(u \in \RR^3, v \in \RR^3\) are the edges of the adjacent vertices
% that defines \(\al\).

W = make_sparse(n,n);
for i=1:3
   i2 = mod(i  ,3)+1;
   i3 = mod(i+1,3)+1;
   u = X(:,F(i2,:)) - X(:,F(i,:));
   v = X(:,F(i3,:)) - X(:,F(i,:));
   % normalize the vectors   
   u = u ./ repmat( sqrt(sum(u.^2,1)), [3 1] );
   v = v ./ repmat( sqrt(sum(v.^2,1)), [3 1] );
   % compute angles
   alpha = 1 ./ tan( acos(sum(u.*v,1)) );
   alpha = max(alpha, 1e-2); % avoid degeneracy
   W = W + make_sparse(F(i2,:),F(i3,:), alpha, n, n );
   W = W + make_sparse(F(i3,:),F(i2,:), alpha, n, n );
end

%%
% Compute the symmetric Laplacian matrix
% \[ L = D-W
%   \qwhereq D = \diag_i\pa{\sum_j W_{i,j}}. \]

d = full( sum(W,1) );
D = spdiags(d(:), 0, n,n);
L = D - W;

%% Fixed Boundary Harmonic Parameterization
% The problem of mesh parameterization corresponds to finding
% 2-D locations \((y_i = (y_i^1,y_i^2) \in \RR^2\) 
% for each original vertex, where \( Y = (y_i)_{i=1}^n \in \RR^{2 \times n} \)
% denotes the flattened positions.

%%
% The goal is for this parameterization to be valid, i.e. the 2-D
% mesh obtained by replacing \(X\) by \(Y\) but keeping the same face connectivity \(F\)
% should not contained flipped faces (all face should have the same orientation in the plane).


%%
% We consider here a linear methods, that finds the parameterization, that
% impose that the coordinates are harmonic inside the domain, and have
% fixed position on the boundary (Dirichlet conditions)
% \[ \forall s=1,2, \quad \forall i \notin B, \quad (L y^s)_i = 0, 
% \qandq \forall j \in B, \quad y^s_j = z_j^s. \]

%%
% In order for this method to define a valid parameterization, it is
% necessary that the fixed position \( z_j = (z^1_j,z^2_j) \in \RR^2 \)
% are consecutive points along a convex polygon.

%%
% Compute the fixed positions \(Z=(z_j)_j\) for the vertices on \(B\).
% Here we use a circle.

p = length(B);
t = linspace(0,2*pi(),p+1); t(p+1) = [];
Z = [cos(t); sin(t)];

%%
% Computing the parameterization requires to solve two independent linear
% system
% \[ \forall s=1,2, \quad  L_1 y^s = r^s \]
% where \(L_1\) is a modified Laplacian, the is obtained from \(L\) by 
% \[ \choice{
%       \forall i \notin B, \quad (L_0)_{i,j} = L_{i,j} \\
%       \forall i \in B, \quad (L_0)_{i,i}=1, \\
%       \forall i \in B, \forall j \neq i, \quad (L_0)_{i,i}=0, 
%   } \]
% i.e. replacing each row indexed by \(B\) by a 1 on the diagonal.

L1 = L;
L1(B,:) = 0;
for i=1:length(B)
    L1(B(i),B(i)) = 1;
end

%%
% Set up the right hand size \(R\) with the fixed position.

R = zeros(2, n);
R(:,B) = Z;

%%
% Solve the two linear systems.
% For Scilab user: this might take some time ...

if using_matlab()
    Y = (L1 \ R')'; 
else
    options.maxit = 300;
    Y(1,:) = perform_cg(L1,R(1,:)',options)';
    Y(2,:) = perform_cg(L1,R(2,:)',options)';
end

%%
% Display the parameterization.

clf;
options.lighting = 0;
plot_mesh([Y;zeros(1,n)],F);
shading('faceted');


%% Mesh Parameterization on a Square
% One can perform a fixed B parameterization on a square.
% This is useful to compute a geometry image (a color image storring the
% position of the vertices).

%EXO
%% Compute the fixed positions \(Z\) of the points indexed by \(B\) that are along a
%% square. Warning: \(p\) is not divisible by 4.
q = floor(p/4);
t = (0:(q-1))/q;
t1 = (0:(p-3*q-1))/(p-3*q);
Z = [[t, t*0+1, 1-t, t1*0]; ...
     [t*0, t, t*0+1, 1-t1]];
clf; hold('on');
hh = plot(Z(1,[1:p 1]),Z(2,[1:p 1]), '.-');
if using_matlab()
    set_linewidth(hh,3);
end
hh = plot(Z(1,[1:p 1]),Z(2,[1:p 1]), 'r.');
if using_matlab()
    set_markersize(hh,20);
end
axis off; axis square; 
%EXO

%EXO
%% Compute the parameterization \(Y\) on a square.
R = zeros(2, n);
R(:,B) = Z;
if using_matlab()
    Y = (L1 \ R')'; 
else
    options.maxit = 300;
    Y(1,:) = perform_cg(L1,R(1,:)',options)';
    Y(2,:) = perform_cg(L1,R(2,:)',options)';
end
clf;
plot_mesh([Y;zeros(1,n)],F);
shading faceted; axis tight;
%EXO

%EXO
%% Shift the \(B\) positions so that the eyes of the model are approximately
%% horizontal.
delta = 6;
sel = [delta:p 1:delta-1];
R = zeros(2,n); R(:,B(sel)) = Z;
% solve
Y = (L1 \ R')'; 
clf;
plot_mesh([Y;zeros(1,n)],F);
shading faceted; axis tight;
%EXO


%% Re-align the Texture
% To map correctly a real image on the surface, the texture needs to be aligned.
% We use here a simple affine mapping to map the eye and mouth of the image
% on the corresponding locaiton on the surface.

%%
% Load a texture image \(T\).

n1 = 256;
T = load_image('lena', n1);
T = T(n1:-1:1,:);

%%
% Display the texture on the mesh, using the parametrization of the mesh as
% texture coordinates.

options.texture = T;
options.texture_coords = Y;
clf;
plot_mesh(X,F,options);
zoom(1.5);

%%
% Position \((u_1,u_2,u_3)\) of the eyes and the mouth in the texture.

u1 = [247;266]*n1/512;
u2 = [247;328]*n1/512;
u3 = [161;301]*n1/512;

%%
% Display.

clf; hold('on');
if using_matlab()
    imageplot(T);
else
    imageplot(T(n1:-1:1,:));
end
h = plot([u1(2) u2(2) u3(2) u1(2)], [u1(1) u2(1) u3(1) u1(1)], '.-b');  
if using_matlab()
    set(h, 'LineWidth', 2);
    set(h, 'MarkerSize', 20);
end

%%
% Positions \((v_1,v_2,v_3)\) of the eyes and the mouth on the parameteric domain

v1 = [310;125]*n1/512;
v2 = [315;350]*n1/512;
v3 = [105;232]*n1/512;

%%
% Display.

clf; hold('on');
plot_graph(W, Y'*(n1-1)+1, 'b');
h = plot([v1(2) v2(2) v3(2) v1(2)], [v1(1) v2(1) v3(1) v1(1)], '.-r');  
axis('off'); axis('equal');
if using_matlab()
    set(h, 'LineWidth', 2);
    set(h, 'MarkerSize', 20);
end


%%
% Apply a similitude to the image so that the mouth and the eye have the correct position in parametric domain.

T1 = perform_image_similitude(T','affine',(u1-1)/n1,(v1-1)/n1,(u2-1)/n1,(v2-1)/n1,(u3-1)/n1,(v3-1)/n1)';

%%
% Display the mesh with the image overlaid.

clf; hold on;
if using_matlab()
    imageplot(T1);
else
    imageplot(T1(n1:-1:1,:));
end
plot_graph(W, Y'*(n1-1)+1, 'b');
h = plot([v1(2) v2(2) v3(2) v1(2)], [v1(1) v2(1) v3(1) v1(1)], '.-r');  
axis('off'); axis('equal');
if using_matlab()
    set(h, 'LineWidth', 2);
    set(h, 'MarkerSize', 20);
end

%%
% Display the texture mapped mesh.

options.texture = T1;
clf;
plot_mesh(X,F,options);
zoom(1.5);

%% Geometry Images
% We use this parameterization to interpolate the mesh on a regular grid.
% This creates a so-called geometry image.

q = 64;
M = zeros(q,q,3);
for i=1:3
    M(:,:,i) = compute_triang_interp(F,Y,X(i,:), q);
end

%%
% Display the re-interpolated mesh as an image.

clf;
imageplot(M);

%%
% Display the re-interpolated mesh as a surface.

C = ones(q,q)*0.8; % color
clf;
colormap(gray(256)); 
surf(M(:,:,1), M(:,:,2), M(:,:,3), C );
view(-40,70);
zoom(1.5);
axis tight; axis square; axis off;
camlight local;  shading faceted;

%%
% Compute a checkboard texture.

[V,U] = meshgrid(1:32,1:32);
T = mod(V+U,2);

%% 
% Display the texture on the surface.

clf;
colormap(gray(256));
plot_surf_texture(M, T);
view(-40,70); zoom(1.5);
axis tight; axis square; axis off;
camlight;


%% Boundary Free Parameterization
% Instead of fixing the positions of the points indexed by \(B\) on a closed convex curve, 
% it is possible to set the boundary free and impose Neumann boundary
% conditions.

%%
% This was proposed in:

%%
% Mathieu Desbrun, Mark Meyer and Pierre Alliez, 
% _Intrinsic Parameterizations of Surface Meshes,_
% In EUROGRAPHICS 2002

%%
% To remove the translation and rotation degree of freedom, one needs to
% impose the position of two vertices.


%%
% Set up a \( (2n,2n) \) linear systeme for the X/Y position of the vertices in
% the plane.

A = make_sparse(2*n,2*n);
A(1:n,1:n) = L;
A(n+1:2*n,n+1:2*n) = L;

%%
% Add the neumann boundary conditions to the system.

for s = 1:p
    s1 = mod(s-2,p)+1;
    s2 = mod(s,p)+1;
    i = B(s);
    i1 = B(s1);
    i2 = B(s2);
    A(i,[i1 i2]+n) = [1, -1];
    A(i+n,[i1 i2]) = [-1, +1];
end

%%
% Select two pinned vertices.

i1 = B(1);
i2 = B(round(p/2));

%%
% Set up the pinned vertices in the system.

A(i1,:) = 0;  A(i1+n,:) = 0;
A(i1,i1) = 1; A(i1+n,i1+n) = 1;
A(i2,:) = 0;  A(i2+n,:) = 0;
A(i2,i2) = 1; A(i2+n,i2+n) = 1;

%%
% Set up the right hand side of the system.

y = zeros(2*n,1);
y([i1 i2]) = [0 1];
y([i1 i2]+n) = [0 0];

%% 
% Solve for the position.

if using_matlab()
    Y = reshape(A\y, [n 2])';
else
    options.maxit = 300;
    Y = reshape( perform_cg(A,y,options), [n 2])';
end
Y = Y - repmat(min(Y,[],2), [1 n]);
Y = Y/max(Y(:));

%%
% Display.

clf;
plot_mesh([Y;zeros(1,n)],F);
shading faceted; axis tight;

%EXO
%% Compute a geometry image from the boundary free parameterization
%% and use it to map a checkboard texture.
q = 64;
M = zeros(q,q,3);
for i=1:3
    M(:,:,i) = compute_triang_interp(F,Y,X(i,:), q);
end
[Y,X] = meshgrid(1:q,1:q);
T = mod(X+Y,2);
clf;
colormap(gray(256));
plot_surf_texture(M, T);
view(-40,70); zoom(1.5);
axis tight; axis square; axis off;
camlight;
%EXO


%% Parameterization of a Larger Mesh
% One can parameterize a larger mesh, having a large area/perimeter of
% B distortion. This usually procudes poor quality parameterization,
% with a large distortion.

%% 
% Load a mesh.

name = 'david-head';
[X,F] = read_mesh(name);
n = size(X,2);
clear options; options.name = name;

%%
% Display the mesh.

clf;
options.lighting = 1;
plot_mesh(X,F,options);
shading faceted;


%EXO
%% Compute the Laplacian \(L\) of the mesh, the boundary \(B\) and the modified
%% Laplacian \(L_1\).
W = make_sparse(n,n);
for i=1:3
   i1 = mod(i-1,3)+1;
   i2 = mod(i  ,3)+1;
   i3 = mod(i+1,3)+1;
   pp = X(:,F(i2,:)) - X(:,F(i1,:));
   qq = X(:,F(i3,:)) - X(:,F(i1,:));
   % normalize the vectors   
   pp = pp ./ repmat( sqrt(sum(pp.^2,1)), [3 1] );
   qq = qq ./ repmat( sqrt(sum(qq.^2,1)), [3 1] );
   % compute angles
   a = 1 ./ tan( acos(sum(pp.*qq,1)) );
   a = max(a, 1e-2); % avoid degeneracy
   W = W + make_sparse(F(i2,:),F(i3,:), a, n, n );
   W = W + make_sparse(F(i3,:),F(i2,:), a, n, n );
end
% Compute the symmetric Laplacian matrix.
d = full( sum(W,1) );
D = spdiags(d(:), 0, n,n);
L = D - W;
%
options.verb = 0;
B = compute_boundary(F, options);
%
L1 = L;
L1(B,:) = 0;
for i=1:length(B)
    L1(B(i),B(i)) = 1;
end
%EXO

%%
% Display the boundary.

clf;
hold on;
plot_mesh(X, F,options);
hh = plot3( X(1,B),X(2,B),X(3,B), 'r' );
set_linewidth(hh,3);
view(-150,-30);

%EXO
%% Perform the parameterization of the mesh on a circle.
p = length(B);
t = linspace(0,2*pi(),p+1)'; t(p) = [];
x0 = cos(t); y0 = sin(t);
% 
Rx = zeros(n,1); Rx(B) = x0;
Ry = zeros(n,1); Ry(B) = y0;
if using_matlab()
    x = L1 \ Rx; 
    y = L1 \ Ry;
else
    options.maxit = 300;
    x = perform_cg(L1,Rx,options);
    y = perform_cg(L1,Ry,options);
end
Y = [x';y'];
clf;
options.lighting = 0;
plot_mesh([Y;zeros(1,n)],F);
shading('faceted');
%EXO


%EXO
%% Display a high frequency function defined on the parameteric domain on
%% the mesh. What do you observe ?
f = 10;
v = cos(f*Y(1,:)*2*pi()) .* cos(f*Y(2,:)*2*pi());
% Display the function on the 2D parameteric domain.
options.face_vertex_color = rescale(v(:));
clf;
subplot(1,2,1);
colormap(jet(256));
plot_mesh([Y;zeros(1,n)],F,options);
colormap(jet(256));
view(2);
subplot(1,2,2);
clf;
colormap(jet(256));
plot_mesh(X,F,options);
colormap(jet(256));
%EXO
