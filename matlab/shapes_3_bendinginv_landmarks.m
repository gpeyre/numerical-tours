%% Geodesic Bending Invariants with Landmarks
% This tour explores the use of farthest point sampling to compute bending
% invariant with classical MDS (strain minimization).

perform_toolbox_installation('signal', 'general', 'graph');

%CMT
rep = 'results/fastmarching_bendinginv_2d/';
if not(exist(rep))
    mkdir(rep);
end
%CMT

%% Farthest Points Landmarks Seeding
% For large mesh, computing all the pairwise distances is intractable.
% It is possible to speed up the computation by restricting the computation
% to a small subset of landmarks.

%%
% This seeding strategy was used for surface remeshing in:

%%
% _Geodesic Remeshing Using Front Propagation_, 
% Gabriel Peyré and Laurent Cohen,
% International Journal on Computer Vision, Vol. 69(1), p.145-156, Aug. 2006.

%%
% Load a mesh.

name = 'elephant-50kv';
options.name = name;
[vertex,faces] = read_mesh(name);
nverts = size(vertex,2);

%%
% Display it.

clf;
plot_mesh(vertex,faces, options);

%CMT
saveas(gcf, [rep name '-original.png'], 'png');
%CMT


%%
% Compute a sparse set of landmarks to speed up the geodesic computations.
% The landmarks are computed using farthest point sampling.

%%
% First landmarks, at random.

landmarks = 23057;
Dland = [];

%%
% Perform Fast Marching to compute the geodesic distance, and record it.

[Dland(:,end+1),S,Q] = perform_fast_marching_mesh(vertex, faces, landmarks(end));

%%
% Select farthest point. Here, |min(Dland,[],2)| is the distance to the set of seed
% points.

[tmp,landmarks(end+1)] = max( min(Dland,[],2) );

%%
% Update distance function.

[Dland(:,end+1),S,Q] = perform_fast_marching_mesh(vertex, faces, landmarks(end));

%%
% Display distances.

clf;
options.start_points = landmarks;
plot_fast_marching_mesh(vertex,faces, min(Dland,[],2) , [], options);

%EXO
%% Compute a set of |n = 300| vertex by iterating this farthest
%% point sampling. Display the progression of the sampling.
n = 400;
landmarks = 1;
Dland = [];
k = 1;
displ = round(linspace(0,1,5)*n); displ(1) = [];
clf;
for i=1:n
    if not(isempty(Dland))
        [tmp,landmarks(end+1)] = max( min(Dland,[],2) );
    end
    [Dland(:,end+1),S,Q] = perform_fast_marching_mesh(vertex, faces, landmarks(end));
    if i==displ(k)
        options.start_points = landmarks;
        subplot(2,2,k);
        options.colorfx = 'equalize';
        plot_fast_marching_mesh(vertex,faces, min(Dland,[],2) , [], options);
        k = k+1;
    end
end
%EXO

%%
% Compute the distance matrix restricted to the landmarks.

D = Dland(landmarks,:);
D = (D+D')/2;


%% Bending Invariant by Strain Minimization and Nistrom Interpolation
% One can compute the bending invariant of the set of landmarks, and then
% apply it to the whole mesh using interpolation.

%%
% Compute a centered kernel for the Landmarks, that should be approximately
% a matrix of inner products.

J = eye(n) - ones(n)/n;
K = -1/2 * J*(D.^2)*J;

%% 
% Perform classical MDS on the reduced set of points, to obtain new positions in 3D.

opt.disp = 0; 
[Xstrain, val] = eigs(K, 3, 'LR', opt);
Xstrain = Xstrain .* repmat(sqrt(diag(val))', [n 1]);
Xstrain = Xstrain';

%%
% Interpolate the locations to the whole mesh by Nystrom
% eigen-extrapolation, as detailed in 

%%
% _Sparse multidimensional scaling using landmark points_
% V. de Silva, J.B. Tenenbaum, Preprint.

vertex1 = zeros(nverts,3);
deltan = mean(Dland.^2,1);
for i=1:nverts
    deltax = Dland(i,:).^2;
    vertex1(i,:) = 1/2 * ( Xstrain * ( deltan-deltax )' )';
end
vertex1 = vertex1';

%%
% Display the bending invariant mesh.

clf;
plot_mesh(vertex1,faces,options);

%CMT
view(50,-65); camlight;
saveas(gcf, [rep name '-strain.png'], 'png');
%CMT       

%% Farthest Point for Stress Minimization
% The proposed interpolation method is valid only for the Strain minimizer
% (spectral Nistrom interpolation). One thus needs to use another
% interpolation method.


%%
% See for instance this work for a method to do such an interpolation:

%%
% A. M. Bronstein, M. M. Bronstein, R. Kimmel, 
% _Efficient computation of isometry-invariant distances between surfaces_, 
% SIAM J. Scientific Computing, Vol. 28/5, pp. 1812-1836, 2006.

%EXO
%% Create an interpolation scheme to interpolate the result of MDS
%% dimensionality reduction with Stree minimization (SMACOF algorithm).
%EXO