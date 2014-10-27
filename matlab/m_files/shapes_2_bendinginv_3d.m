%% Geodesic Bending Invariants for Surfaces
% This tour explores the computation of bending invariants of surfaces.

perform_toolbox_installation('signal', 'general', 'graph');

%CMT
rep = 'results/fastmarching_bendinginv_3d/';
if not(exist(rep))
    mkdir(rep);
end
%CMT

%% Bending Invariants
% Bending invariants replace the position of the vertices in a shape \(\Ss\) (2-D or 3-D)
% by new positions that are insensitive to isometric deformation of the shape.
% This defines a bending invariant signature that can be used for surface
% matching.
  
%% 
% Bending invariant were introduced in <#biblio [EladKim03]>.
% A related method was developped for brain flattening in <#biblio [SchwShWolf89]>.
% This method is related to the Isomap algorithm for manifold learning
% <#biblio [TenSolvLang03]>.

%% 
% We assume that \(Ss\) has some Riemannian metric, for instance coming
% from the embedding of a surface in 3-D Euclidian space, or by restriction of 
% the Euclian 2-D space to a 2-D sub-domain (planar shape). One thus can
% compute the geodesic distance \(d(x,x')\) between points \(x,x' \in
% \Ss\).

%%
% The bending invariant \(\tilde \Ss\) of \(\Ss\) is defined as the set of
% points \(Y = (y_i)_j \subset \RR^d\) that are optimized so that the Euclidean
% distance between points in \(Y\) matches as closely the geodesic distance
% between points in \(X\), i.e.
% \[ \forall i, j, \quad \norm{y_i-y_j} \approx d(x_i,x_j) \]

%%
% Multi-dimensional scaling (MDS) is a class of method that aims at
% computing such a set of points \(Y \in \RR^{d \times N}\) in \(\RR^d\)
% such that 
% \[ \forall i, j, \quad \norm{y_i-y_j} \approx \de_{i,j} \]
% where \(\de \in \RR^{N \times N}\) is a input data matrix. 
% For a detailed overview of MDS, we refer to the book <#biblio [BorgGroe97]>

%%
% In this tour, we apply two specific MDS algorithms (strain and stress
% minimization) with input \(\de_{i,j} = d(x_i,x_j)\).

%% 3-D Surfaces and Geodesic Distances
% We consider here a syrface \(\Ss \subset \RR^3\).

%%
% Load a mesh of \(N\) vertices that discretizes this surfaces.

name = 'camel';
options.name = name;
[V,F] = read_mesh(name);
N = size(V,2);

%%
% Display it.

clf;
plot_mesh(V,F, options);

%CMT
saveas(gcf, [rep name '-original.png'], 'png');
%CMT

%%
% The geodesic distance map \(U(x) = d(x,x_i)\) to a starting point \(x_i\)
% can be computed in \(O(N \log(N))\) operations on a mesh of \(N\) 
% vertices using the Fast Marching algorithm.

i = 1;
U = perform_fast_marching_mesh(V, F, i);

%%
% Extract a bunch of geodesic shortest paths from \(x_i\) 
% to randomly selected vertices \( (x_j)_{j \in J} \).

options.method = 'continuous';
J = randperm(N); J = J(1:50);
paths = compute_geodesic_mesh(U, V, F, J, options);

%% 
% Display the distance \(U\) on the 3-D mesh
% together with the geodesic paths.

clf;
plot_fast_marching_mesh(V, F, U, paths, options);

%%
% The geodesic distance matrix \(\de \in \RR^{N \times N}\) is defined as
% \[ \forall i,j=1,\ldots,N, \quad
%       \de_{i,j} = d(x_i,x_j). \]

%EXO
%% Compute the geodesic distance matrix \(\de\).
%% It is going to take some of time.
delta = zeros(N,N);
for i=1:N
%    progressbar(i,N);
    [delta(:,i),S,Q] = perform_fast_marching_mesh(V, F, i);
end
delta = (delta+delta)'/2;
%EXO


%% Bending Invariant with Strain Minimization
% The goal is to compute a set of points \(Y = (y_i)_{i=1}^N\) in
% \(\RR^d\), (here we use \(d=2\)) stored in a matrix \(Y \in \RR^{d \times N}\)
% such that
% \[ \forall i, j, \quad D^2(Y)_{i,j} \approx \de_{i,j}^2 
%   \qwhereq  D^2(Y)_{i,j} = \norm{y_i-y_j}^2. \]

%%
% Target dimensionality \(d\).

d = 3;

%%
% This can be achieved by minimzing a \(L^2\) loss
% \[ \umin{Y} \norm{ D^2(Y)-\de^2 }^2 = 
%       \sum_{i<j} \abs{ \norm{y_i-y_j}^2 - \de_{i,j}^2 }^2. \]

%%
% Strain minimization consider instead the following weighted \(L^2\) loss
% (so-called strain)
% \[ \umin{Y \in \RR^{d \times N} } \text{Strain}(Y) 
%       = \norm{ J ( D^2(Y)-\de^2 ) J }^2
% \]
% where \(J\) is the so-called centering matrix
% \[ J_{i,j} = 
%       \choice{
%           1-1/N \qifq i=j, \\
%           -1/N \qifq i \neq j.
%  }\]

J = eye(N) - ones(N)/N;

%%
% Using the properties of squared-distance matrices \(D^2(Y)\), one can
% show that 
% \[ \norm{ J ( D^2(Y)-\de^2 ) J }^2 = 
%   \norm{ Y Y^* - K }^2 
%   \qwhereq K = - \frac{1}{2} J \de^2 J. \]

K = -1/2 * J*(delta.^2)*J;

%%
% The solution to this (non-convex) optimization problem can be computed
% exactly as the rows of \(Y\) being the two leading eigenvectors of \(K\)
% propery rescaled.

opt.disp = 0; 
[Y, v] = eigs(K, d, 'LR', opt);
Y = Y .* repmat(sqrt(diag(v))', [N 1]);
Y = Y';

%%
% Display the bending invariant surface.

clf;
plot_mesh(Y,F, options);

%CMT
saveas(gcf, [rep name '-bendinv-strain.png'], 'png');
%CMT


%% Bending Invariant with Stress Minimization
% The stress functional does not have geometrical meaning. 
% An alternative MDS method directly minimizes a geometric loss, the
% so-called Stress
% \[ \umin{Y \in \RR^{d \times N} } \text{Stress}(Y) = 
%       \norm{ D(Y)-\de }^2 = 
%       \sum_{i<j} \abs{ \norm{y_i-y_j} - \de_{i,j} }^2. 
%  \]
% It is possible to find a local minimizer of this energy by various
% descent algorithms, as initially proposed by <#biblio [Kruskal64]>

Stress = @(d)sqrt( sum( abs(delta(:)-d(:)).^2 ) / N^2 ); 

%%
% Operator to compute the distance matrix \(D(Y)\).

D = @(Y)sqrt( repmat(sum(Y.^2),N,1) + repmat(sum(Y.^2),N,1)' - 2*Y'*Y);

%%
% The SMACOF (Scaling by majorizing a convex function) algorithm 
% solves at each iterations a quadratic energy, that is guaranteed to 
% diminish the value of the Strain. It was introduced by <#biblio [Leeuw77]>

%%
% It computes iterates \(X^{(\ell)}\) as
% \[ X^{(\ell+1)} = \frac{1}{N} X^{(\ell)} B(D(X^{(\ell)}))^*,  \]
% where 
% \[ B(D) = \choice{
%          -\frac{\de_{i,j}}{D_{i,j}} \qifq i \neq j, \\
%          -\sum_{k} B(D)_{i,k} \qifq i = j.
% } \]

%%
% Initialize the positions for the algorithm.

Y = V;

%%
% Operator \(B\).

remove_diag = @(b)b - diag(sum(b));
B = @(D1)remove_diag( -delta./max(D1,1e-10) );

%%
% Update the positions.

Y = Y * B(D(Y))' / N;

%EXO
%% Perform the SMACOF iterative algorithm.
%% Save in a variable |s(l)| the values of
%% Stress\(( X^{(\ell)} )\).
niter = 50;
stress = [];
Y = V;
ndisp = [1 5 10 niter Inf];
s = [];
k = 1;
clf;
for i=1:niter
    if ndisp(k)==i
        subplot(2,2,k); 
        plot_mesh(Y,F, options);
        axis('equal'); axis('off');
        k = k+1;
    end
    Y = Y * B(D(Y))' / N;  
    % update
    % Y = Y-repmat(mean(Y,2), [1 N]);
    % record stress
    s(end+1) = Stress(D(Y));
end
axis('equal'); axis('off'); axis('ij');
%EXO

%%
% Plot stress evolution during minimization.

clf;
plot(s, '.-', 'LineWidth', 2, 'MarkerSize', 20);
axis('tight');

%%
% Plot the optimized invariant shape.

clf;        
plot_mesh(Y,F, options);

%% Surface Retrieval with Bending Invariant.
% One can compute a bending invariant signature for each mesh in a library
% of 3D surface.

%%
% Isometry-invariant retrival is then perform by matching the invariant
% signature.

%EXO
%% Implement a surface retrival algorithm based on these bending invariants.
% No correction for this exercise.
%EXO



%% Bibliography
% <html><a name="biblio"></a></html>

%%
% * [EladKim03] A. Elad and R. Kimmel, <http://dx.doi.org/10.1109/TPAMI.2003.1233902 _On bending invariant signatures for surfaces_>, IEEE Transactions onPattern Analysis and Machine Intelligence, Vol. 25(10), p. 1285-1295, 2003.
% * [SchwShWolf89] E.L. Schwartz and A. Shaw and E. Wolfson, <http://dx.doi.org/10.1109/34.35506 _A Numerical Solution to the Generalized Mapmaker's Problem: Flattening Nonconvex Polyhedral Surfaces_>, IEEE Transactions on Pattern Analysis and Machine Intelligence, 11(9), p. 1005-1008, 1989.
% * [TenSolvLang03] J. B. Tenenbaum, V. de Silva and J. C. Langford, <http://dx.doi.org/10.1126/science.290.5500.2319 _A Global Geometric Framework for Nonlinear Dimensionality Reduction_>, Science 290 (5500): 2319-2323, 22 December 2000 
% * [Kruskal64] J. B. Kruskal, <http://dx.doi.org/10.1007/BF02289565 _Multidimensional scaling by optimizing goodness of fit to a nonmetric hypothesis_>, Psychometrika 29 (1): 1?27, 1964.
% * [Leeuw77] J. de Leeuw, <http://statistics.ucla.edu/preprints/uclastat-preprint-2005:35 _Applications of convex analysis to multidimensional scaling_>, in Recent developments in statistics, pp. 133?145, 1977
% * [BorgGroe97] I. Borg and P. Groenen, <http://www.springeronline.com/0-387-25150-2 _Modern Multidimensional Scaling: theory and applications_>, New York: Springer-Verlag, 1997.

