%% Heuristically Driven Front Propagation
% This tour explores the use of heuristics to speed up Fast Marching methods in 2D.

%%
% The use of heuristics for the computation of geodesic paths was
% introduced in 

%% 
% _Heuristically Driven Front Propagation for Fast Geodesic Extraction_
% Gabriel Peyre and Laurent Cohen, 
% International Journal for Computational Vision and Biomechanics, Vol. 1(1), p.55-67, Jan-June 2008.

perform_toolbox_installation('signal', 'general', 'graph');

%CMT
rep = 'results/fastmarching_heuristics/';
if not(exist(rep))
    mkdir(rep);
end
%CMT


%% Ideal Heuristically Driven Front Propagation
% The ideal heuristic |H| is the remaining distance to the ending point.
% It is ideal in the sense that it is as difficult to compute this distance
% than to solve for the original problem of extracting the geodesic.

%%
% One can however study the influence of this heuristic by replacing |H|
% by |weight*H| where |weight<1| is a sub-optimality factor.

%%
% Load an image.

n = 300;
name = 'road2';
M = rescale(load_image(name,n));

%%
% Display it.

clf;
imageplot(M);

%%
% Define start and end points (you can use your own points).

pstart = [14;161];
pend = [293;148];

%%
% Compute a metric to extact the road.

W = abs(M-M(pstart(1),pstart(2)));
W = rescale(W, 1e-2,1);

%%
% Display it.

clf;
imageplot(W);


%%
% Perform the full propagation.

[D,S] = perform_fast_marching(1./W, pstart);

%% 
% Extract a geodesic curve.

p = compute_geodesic(D,pend);

%%
% Display the distance and the geodesic curve.

clf; hold on;
imageplot(convert_distance_color(D,M), 'Distance');
h = plot(p(2,:),p(1,:), '.k'); set(h, 'LineWidth', 2);
h = plot(pstart(2),pstart(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(pend(2),pend(1), '.b'); set(h, 'MarkerSize', 25);
axis ij;

%%
% Compute the heuristic.

[H,S] = perform_fast_marching(1./W, pend);

%%
% Display the ideal heuristic function.

clf; hold on;
imageplot(convert_distance_color(H,M), 'Distance');
h = plot(p(2,:),p(1,:), '.k'); set(h, 'LineWidth', 2);
h = plot(pstart(2),pstart(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(pend(2),pend(1), '.b'); set(h, 'MarkerSize', 25);
axis ij;

%EXO
%% Display the set of points satisfying |D+H<=T| for several value of the
%% threshold |T>=D(pend)|. What do you observe ?
d = (D(pend(1),pend(2)) + H(pstart(1),pstart(2)))/2;
Tlist = [1.01 1.05 1.1 1.2] * d;
clf;
for i=1:4
    t = Tlist(i);
    U = cat(3,M,M,M);
    I = find( D+H<=t );
    U(I) = 1; U([I+n^2, I+2*n^2]) = U([I+n^2, I+2*n^2])*.3;
    subplot(2,2,i);
    hold on;
    imageplot(U);
    h = plot(p(2,:),p(1,:), '.k'); set(h, 'LineWidth', 2);
    h = plot(pstart(2),pstart(1), '.g'); set(h, 'MarkerSize', 25);
    h = plot(pend(2),pend(1), '.b'); set(h, 'MarkerSize', 25);
    axis ij;
end
%EXO

%CMT
for i=1:4
    t = Tlist(i);
    U = cat(3,M,M,M);
    I = find( D+H<=t );
    U(I) = 1; U([I+n^2, I+2*n^2]) = U([I+n^2, I+2*n^2])*.3;
    clf;
    hold on;
    imageplot(U);
    h = plot(p(2,:),p(1,:), '.k'); set(h, 'LineWidth', 2);
    h = plot(pstart(2),pstart(1), '.g'); set(h, 'MarkerSize', 25);
    h = plot(pend(2),pend(1), '.b'); set(h, 'MarkerSize', 25);
    axis ij;
    saveas(gcf, [rep name '-ellipsoid-' num2str(i) '.eps'], 'epsc');
end
%CMT

%%
% Perform the heuristically driven propagation.

weight = .9;
options.end_points = pend;
options.heuristic = weight*H;
options.nb_iter_max = Inf;
options.constraint_map = Inf+zeros(n);
[D,S] = perform_fast_marching(1./W, pstart, options);

%%
% Display the region explored by the algorithm.

I = find(S<0);
U = cat(3,M,M,M);
U(I) = 1; U([I+n^2, I+2*n^2]) = U([I+n^2, I+2*n^2])*.3;
clf; hold on;
imageplot(U);
h = plot(p(2,:),p(1,:), '.k'); set(h, 'LineWidth', 2);
h = plot(pstart(2),pstart(1), '.g'); set(h, 'MarkerSize', 25);
h = plot(pend(2),pend(1), '.b'); set(h, 'MarkerSize', 25);
axis ij;

%EXO
%% Display the explored region for different values of |weight|.
wlist = [.6 .8 .9 .95];
clf;
for i=1:4
    weight = wlist(i);
    options.end_points = pend;
    options.heuristic = weight*H;
    options.nb_iter_max = Inf;
    options.constraint_map = Inf+zeros(n);
    [D,S] = perform_fast_marching(1./W, pstart, options);
    %
    I = find(S<0);
    U = cat(3,M,M,M);
    U(I) = 1; U([I+n^2, I+2*n^2]) = U([I+n^2, I+2*n^2])*.3;
    subplot(2,2,i); 
    hold on;
    imageplot(U);
    h = plot(p(2,:),p(1,:), '.k'); set(h, 'LineWidth', 2);
    h = plot(pstart(2),pstart(1), '.g'); set(h, 'MarkerSize', 25);
    h = plot(pend(2),pend(1), '.b'); set(h, 'MarkerSize', 25);
    axis ij;
end
%EXO

%CMT
for i=1:4
    weight = wlist(i);
    options.end_points = pend;
    options.heuristic = weight*H;
    options.nb_iter_max = Inf;
    options.constraint_map = Inf+zeros(n);
    [D,S] = perform_fast_marching(1./W, pstart, options);
    %
    I = find(S<0);
    U = cat(3,M,M,M);
    U(I) = 1; U([I+n^2, I+2*n^2]) = U([I+n^2, I+2*n^2])*.3;
    clf; hold on;
    imageplot(U);
    h = plot(p(2,:),p(1,:), '.k'); set(h, 'LineWidth', 2);
    h = plot(pstart(2),pstart(1), '.g'); set(h, 'MarkerSize', 25);
    h = plot(pend(2),pend(1), '.b'); set(h, 'MarkerSize', 25);
    axis ij;
    saveas(gcf, [rep name '-ideal-heuristic-' num2str(i) '.eps'], 'epsc');
end
%CMT



%% Landmark-based Heuristically Driven Front Propagation
% An heuristic can be derived using a pre-computed set of distance to
% landmark points. The more landmark, the more accurate the heuristic is.

%%
% Compute randomized landmarks.

q = 10;
landmarks = floor(rand(2,q)*n)+1;

%%
% Pre-compute distances to landmarks.

Dland = zeros(n,n,q);
for i=1:q
    Dland(:,:,i) = perform_fast_marching(1./W, landmarks(:,i));
end

%%
% Compute the heuristic.

Dend = Dland( pend(1), pend(2), :);
H = max(abs(Dland-repmat(Dend, [n n 1])), [], 3);

%%
% Display the heuristic.

clf;
hold on;
imageplot(H);
contour(H, 10, 'k', 'LineWidth', 2);
colormap jet(256);
h = plot(landmarks(1,:), landmarks(2,:), 'y.');
set(h, 'MarkerSize', 15);
axis ij;

%%
% That should be compared with the optimal ideal heuristic.

[H0,S] = perform_fast_marching(1./W, pend);
clf;
hold on;
imageplot(H0);
contour(H0, 10, 'k', 'LineWidth', 2);
colormap jet(256);
axis ij;

%EXO
%% Display the convergence of the heuristic as the number of landmark
%% increases.
qlist = [2 5 10 50];
q = max(qlist);
landmarks = floor(rand(2,q)*n)+1;
Dland = zeros(n,n,q);
for i=1:q
    Dland(:,:,i) = perform_fast_marching(1./W, landmarks(:,i));
end
clf;
for i=1:4
    q = qlist(i);
    Dend = Dland( pend(1), pend(2), :);
    H = max(abs(Dland(:,:,1:q)-repmat(Dend(1:q), [n n 1])), [], 3);
    subplot(2,2,i);
    hold on;
    imageplot(H);
    contour(H, 10, 'k', 'LineWidth', 2);
    colormap jet(256);
    h = plot(landmarks(1,1:q), landmarks(2,1:q), 'y.');
    set(h, 'MarkerSize', 15);
    axis ij;
end
%EXO

%CMT
clf;
for i=1:4
    q = qlist(i);
    Dend = Dland( pend(1), pend(2), :);
    H = max(abs(Dland(:,:,1:q)-repmat(Dend(1:q), [n n 1])), [], 3);
    clf;
    hold on;
    imageplot(H);
    contour(H, 10, 'k', 'LineWidth', 2);
    colormap jet(256);
    h = plot(landmarks(1,1:q), landmarks(2,1:q), 'y.');
    set(h, 'MarkerSize', 15);
    axis ij;
    saveas(gcf, [rep name '-landmark-heuristic-' num2str(i) '.eps'], 'epsc');
end
%CMT

%EXO
%% Perform the heuristically driven propagation with a landmark-based
%% heuristic.
clf;
for i=1:4
    q = qlist(i);
    Dend = Dland( pend(1), pend(2), :);
    H = max(abs(Dland(:,:,1:q)-repmat(Dend(1:q), [n n 1])), [], 3);
    %
    options.end_points = pend;
    options.heuristic = H;
    options.nb_iter_max = Inf;
    options.constraint_map = Inf+zeros(n);
    [D,S] = perform_fast_marching(1./W, pstart, options);
    %
    I = find(S<0);
    U = cat(3,M,M,M);
    U(I) = 1; U([I+n^2, I+2*n^2]) = U([I+n^2, I+2*n^2])*.3;
    subplot(2,2,i); 
    hold on;
    imageplot(U);
    h = plot(p(2,:),p(1,:), '.k'); set(h, 'LineWidth', 2);
    h = plot(pstart(2),pstart(1), '.g'); set(h, 'MarkerSize', 25);
    h = plot(pend(2),pend(1), '.b'); set(h, 'MarkerSize', 25);
    h = plot(landmarks(1,1:q), landmarks(2,1:q), 'y.'); set(h, 'MarkerSize', 15);
    axis ij;
end
%EXO

%CMT
for i=1:4
    q = qlist(i);
    Dend = Dland( pend(1), pend(2), :);
    H = max(abs(Dland(:,:,1:q)-repmat(Dend(1:q), [n n 1])), [], 3);
    %
    options.end_points = pend;
    options.heuristic = H;
    options.nb_iter_max = Inf;
    options.constraint_map = Inf+zeros(n);
    [D,S] = perform_fast_marching(1./W, pstart, options);
    %
    I = find(S<0);
    U = cat(3,M,M,M);
    U(I) = 1; U([I+n^2, I+2*n^2]) = U([I+n^2, I+2*n^2])*.3;
    clf; 
    hold on;
    imageplot(U);
    h = plot(p(2,:),p(1,:), '.k'); set(h, 'LineWidth', 2);
    h = plot(pstart(2),pstart(1), '.g'); set(h, 'MarkerSize', 25);
    h = plot(pend(2),pend(1), '.b'); set(h, 'MarkerSize', 25);
    h = plot(landmarks(1,1:q), landmarks(2,1:q), 'y.'); set(h, 'MarkerSize', 15);
    axis ij;
    saveas(gcf, [rep name '-landmark-extraction-' num2str(i) '.eps'], 'epsc');
end
%CMT

%EXO
%% Find a strategy to find optimal seeding position for the landmarks.
%EXO


%% Heuristics on 3D Meshes
% It is possible to use the same heuristics to drive the computation of
% geodesic paths on 3D meshes.

%EXO
%% Perform the landmark-based heuristically driven propagation on a 3D
%% mesh.
%EXO
