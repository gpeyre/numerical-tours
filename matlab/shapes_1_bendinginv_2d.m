%% Geodesic Bending Invariants for Shapes
% This tour explores the computation of bending invariants of shapes.

perform_toolbox_installation('signal', 'general', 'graph');

%CMT
rep = 'results/fastmarching_bendinginv_2d/';
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

%% 2-D Shapes
% We consider here the case where \(\Ss\) is a sub-domain of \(\RR^2\).

%%
% A binary shape \(\Ss\) is represented as a binary image \(S \in \RR^Q\)
% of \(Q=q \times q\) pixels.

clear options;
q = 400;
name = 'centaur1';
S = load_image(name,q);
S = perform_blurring(S,5);
S = double( rescale( S )>.5 );
if S(1)==1 
    S = 1-S;
end

%%
% Compute its boundary \(b = (b_i)_{i=1}^L \in \RR^{2 \times L}\).

b = compute_shape_boundary(S);
L = size(b,2);

%%
% Display the shape.

lw = 2;
clf; hold on;
imageplot(-S);
plot(b(2,:), b(1,:), 'r', 'LineWidth', lw); axis ij;


%% Geodesic Distance
% We consider the geodesic distance obtained by contraining the shortest
% curve to be inside \(\Ss\)
% \[ d(x,x') = \umin{ \ga } \int_0^1 \abs{\ga'(t)} d t, \]
% where \(\ga\) should satisfy
% \[ \ga(0)=x, \quad \ga(1)=x' \qandq \forall t \in [0,1], \: \ga(t) \in \Ss. \]

%%
% The geodesic distance map \(U(x) = d(x,x_0)\) to a starting point \(x_0\)
% can be computed in \(O(q^2 \log(q))\) operations on a grid of \(q \times
% q\) points using the Fast Marching algorithm.


%%
% Design a constraint map for the Fast-Marching, to enforce the propagation
% inside the shape.

S1 = perform_convolution(S, ones(3)/9)<.01; S1=double(S1==0);
C = zeros(q)-Inf; C(S1==1) = +Inf;
options.constraint_map = C;

%% 
% Shortcut to compute the geodesic distance to some point \(x\).

geod = @(x)perform_fast_marching(ones(q), x, options);

%%
% Compute the geodesic distance to a point \(x\)

x = [263; 55];
options.nb_iter_max = Inf;
U = geod(x);

%%
% Display the distance and geodesic curves. 

clf;
options.display_levelsets = 1;
options.pstart = x;
options.nbr_levelsets = 60;
U(S==0) = Inf;
display_shape_function(U, options);

%EXO
%% Compute curves joining the start point to several points along the
%% boundary.
ms = 20;
U = geod(x);
npaths = 30;
sel = round(linspace(1,L+1,npaths+1)); sel(end) = [];
end_points = b(:,sel);
clf; hold on;
imageplot(1-S);
for i=1:npaths
    p = compute_geodesic(U,end_points(:,i));
    h = plot( p(2,:), p(1,:), 'g' ); set(h, 'LineWidth', lw);
end
h = plot(x(2),x(1), '.r'); set(h, 'MarkerSize', ms);
h = plot(end_points(2,:),end_points(1,:), '.b'); set(h, 'MarkerSize', ms);
axis ij;
%EXO

%% Geodesic Distance Matrix
% We define a set \(X = (x_i)_{i=1}^n \in \RR^{2 \times N} \subset \Ss\)
% of sampling points.

%%
% Total number \(N\) of sampling points.

N = 1000;

%%
% The sampling is made of \(N_0\) points on the boundary \(b\), 
% and \(N-N_0\) points inside the shape \(\Ss\).

%%
% Number of points on the boundary. 

N0 = round(.4*N);

%%
% Sampling on the boundary. 

I = round(linspace(1,L+1,N0+1));
X = round(b(:,I(1:end-1)));

%% 
% Add \(N-N_0\) points inside the shape.

[y,x] = meshgrid(1:q,1:q);
I = find(S==1);
I = I(randperm(length(I))); I = I(1:N-N0);
X(:,end+1:N) = [x(I),y(I)]';

%%
% Display the sampling points.

clf; hold on;
imageplot(1-S);
plot(X(2,:), X(1,:), 'r.', 'MarkerSize', 15);
axis('ij');

%%
% The geodesic distance matrix \(\de \in \RR^{N \times N}\) is defined as
% \[ \forall i,j=1,\ldots,N, \quad
%       \de_{i,j} = d(x_i,x_j). \]

%EXO
%% Compute the geodesic distance matrix \(\de\).
delta = zeros(N,N);
sel = X(1,:) + (X(2,:)-1)*q;
for i=1:N
    U = geod(X(:,i));
    delta(:,i) = U(sel);
end
delta = (delta+delta')/2;
%EXO


%%
% Display the geodesic matrix for \(1 \leq i,j \leq N_0\)
% (points along the boundary \(b\)).

clf;
imageplot(delta(1:N0,1:N0));
colormap(jet(256));

%% Bending Invariant with Strain Minimization
% The goal is to compute a set of points \(Y = (y_i)_{i=1}^N\) in
% \(\RR^d\), (here we use \(d=2\)) stored in a matrix \(Y \in \RR^{d \times N}\)
% such that
% \[ \forall i, j, \quad D^2(Y)_{i,j} \approx \de_{i,j}^2 
%   \qwhereq  D^2(Y)_{i,j} = \norm{y_i-y_j}^2. \]

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
[Y, v] = eigs(K, 2, 'LR', opt);
Y = Y .* repmat(sqrt(diag(v))', [N 1]);
Y = Y';

%%
% Extract the boundary part of the mapped data.
% Rotate the shape if necessary.

theta = -.8*pi;
uv = Y(:,1:N0);
uv = [cos(theta)*uv(1,:) + sin(theta)*uv(2,:); - sin(theta)*uv(1,:) + cos(theta)*uv(2,:)];

%%
% Display the bending invariant boundary.

clf;
h = plot(uv(2,:), uv(1,:)); 
axis('ij'); axis('equal'); axis('off');
set(h, 'LineWidth', lw);

%CMT
clf;
h = plot(X(2,[1:N0 1]), X(1,[1:N0 1]), 'k'); 
axis('equal'); axis('off'); axis('ij');
set(h, 'LineWidth', 2);
saveas(gcf, [rep 'shape-original.eps'], 'eps');
d = sqrt(sum(uv.^2)); I = find(d<=1.5*mean(d));
clf;
h = plot(uv(2,I), uv(1,I), 'k'); 
axis('equal'); axis('off'); axis('ij');
set(h, 'LineWidth', 2);
saveas(gcf, [rep 'shape-strain.eps'], 'eps');
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

Y = X/q;

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
s = [];
Y = X/q;
ndisp = [1 5 niter Inf];
clf; k =1; 
hold on;
Y = Y-repmat(mean(Y,2), [1 N]);
h = plot(Y(2,[1:N0 1]), Y(1,[1:N0 1]));
for i=1:niter
    Y = Y * B(D(Y))' / N;  
    % update
    Y = Y-repmat(mean(Y,2), [1 N]);
    % record stress
    s(end+1) = Stress(D(Y));
    if ndisp(k)==i
        plot(Y(2,[1:N0 1]), Y(1,[1:N0 1])); 
        axis('equal'); axis('off');
        k = k+1;
    end
end
axis('equal'); axis('off'); axis('ij');
%EXO

%CMT
niter = 50;
stress = [];
Y = X/q;
ndisp = [1 5 10 niter Inf];
k =1; 
Y = Y-repmat(mean(Y,2), [1 N]);
clf;
h = plot(Y(2,[1:N0 1]), Y(1,[1:N0 1]), 'k');
axis('equal'); axis('off'); axis('ij');
set(h, 'LineWidth', 2);
saveas(gcf, [rep 'shape-smacof-' num2str(0) '.eps'], 'eps');
for i=1:niter
    % Compute the distance matrix.
    D1 = repmat(sum(Y.^2,1),N,1);
    D1 = sqrt(D1 + D1' - 2*Y'*Y);
    % Compute the scaling matrix.
    B = -delta./max(D1,1e-10);
    B = B - diag(sum(B));
    % update
    Y = (J*(B*Y'))' / N;
    Y = Y-repmat(mean(Y,2), [1 N]);
    % record stress
    stress(end+1) = sqrt( sum( abs(delta(:)-D1(:)).^2 ) / N^2 );
    if ndisp(k)==i
        clf;
        h = plot(Y(2,[1:N0 1]), Y(1,[1:N0 1]), 'k'); 
        axis('equal'); axis('off'); axis('ij');
        set(h, 'LineWidth', 2);
        saveas(gcf, [rep 'shape-smacof-' num2str(k) '.eps'], 'eps');
        k = k+1;
    end
end
%CMT

%%
% Plot stress evolution during minimization.

clf;
plot(s, '.-', 'LineWidth', 2, 'MarkerSize', 20);
axis('tight');

%%
% Plot the optimized invariant shape.

clf;
h = plot(Y(2,[1:N0 1]), Y(1,[1:N0 1])); 
axis('ij'); axis('equal'); axis('off');
set(h, 'LineWidth', lw);


%% Shape Retrieval with Bending Invariant.
% One can compute a bending invariant signature for each shape in a library
% of images.

%%
% Isometry-invariant retrival is then perform by matching the invariant
% signature.

%EXO
%% Implement a shape retrival algorithm based on these bending invariants.
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
