%% Geodesic Farthest Point Sampling
% This tour explores the use geodesic computation to perform image
% sampling.

perform_toolbox_installation('signal', 'general', 'graph');

%CMT
rep = 'results/fastmarching_sampling/';
if not(exist(rep))
    mkdir(rep);
end
%CMT

%% Voronoi Segmentation
% A geodesic Voronoi segementation is obtained by computing a geodesic distance from
% multiple starting points.

%%
% Compute an image with bumps.

n = 256; % size of the image
sigma = n/8; % width of the bumps
[Y,X] = meshgrid(1:n,1:n);
x = n/4; y = n/4;
M = exp( -( (X-x).^2 + (Y-y).^2 )/(2*sigma^2) );
x = 3*n/4; y = 3*n/4;
M = M + exp( -( (X-x).^2 + (Y-y).^2 )/(2*sigma^2) );

%%
% Compute a metric by rescaling |M|.

W = rescale(M,1e-2,1);

%%
% Create random starting points.

m = 20; 
pstart = floor( rand(2,m)*(n-1) ) +1;


%%
% Perform the propagation using the Fast Marching.

[D,Z,Q] = perform_fast_marching(1./W, pstart);


%%
% Display the distance map.

clf;
subplot(1,2,1);
hold on;
imageplot(perform_hist_eq(D,'linear'));  title('Geodesic distance'); 
plot(pstart(2,:), pstart(1,:), 'r.');
subplot(1,2,2);
hold on;
imageplot(Q); title('Voronoi'); 
plot(pstart(2,:), pstart(1,:), 'r.');
colormap jet(256);

%% Geodesic Delaunay Triangulation
% A geodesic Delaunay triangulation is obtained by linking starting points
% whose Voronoi cells touch. This is the dual of the original Voronoi segmentation.

%EXO
%% Using |Q|, compute the faces |faces| of the Delaunay triangulation. 
%% To that end, 
%% extract each quad of values |Q(i,j),Q(i+1,j),Q(i+1,j+1),Q(i,j+1)|, 
%% and add a new face when three of these four values are different (this
%% corresponds to a Voronoi point).
%% Display the obtained triangulation.
% compute quad of values
V = [];
v = Q(1:end-1,1:end-1); V = [V v(:)];
v = Q(2:end,1:end-1); V = [V v(:)];
v = Q(1:end-1,2:end); V = [V v(:)];
v = Q(2:end,2:end); V = [V v(:)];
V = sort(V,2);
V = unique(V, 'rows');
d = (V(:,1)~=V(:,2)) + (V(:,2)~=V(:,3)) + (V(:,3)~=V(:,4));
V = V(d==2,:);
for i=1:size(V,1)
    V(i,1:3) = unique(V(i,:));
end
faces = V(:,1:3)';
%EXO

%%
% Display the obtained triangulation.

clf;
subplot(1,2,1);
hold on;
imageplot(Q, 'Voronoi'); axis ij;
plot(pstart(2,:), pstart(1,:), 'r.');
subplot(1,2,2);
plot_triangulation(pstart,faces, M);
colormap jet(256);

%%
% Compare this Geodesic Delaunay triangulation with the Euclidean
% triangulation.

faces_euc = compute_delaunay(pstart);

%%
% Display.

clf;
subplot(1,2,1);
plot_triangulation(pstart,faces, M); title('Geodesic');
subplot(1,2,2);
plot_triangulation(pstart,faces_euc, M); title('Euclidean');
colormap jet(256);

%% Farthest Point Sampling
% To sample point uniformly according to the geodesic distance, one can use
% an iterative farthest point sampling scheme.

%% 
% Construct a metric that is large in area where we want to sample more
% points. The front should move slowly in high density sampling region.

W = rescale(M,3*1e-1,1);

%%
% Choose the initial point.

vertex = [1;1];

%%
% Compute the geodesic distance.

[D,Z,Q] = perform_fast_marching(1./W, vertex);

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
imageplot(W, 'Metric'); axis ij;
plot(vertex(2,1), vertex(1,1), 'r.');
plot(vertex(2,2), vertex(1,2), 'b.');
subplot(1,2,2);
hold on;
imageplot( perform_hist_eq(D, 'linear'), 'Distance'); axis ij;
plot(vertex(2,1), vertex(1,1), 'r.');
plot(vertex(2,2), vertex(1,2), 'b.');
colormap jet(256);

%%
% Update the value of the distance map with a partial propagation from the last added point.

options.constraint_map = D;
[D1,Z,Q] = perform_fast_marching(1./W, vertex(:,end), options);
% display old/new
clf;
imageplot( convert_distance_color(D,M,1), 'Update distance', 1,3,1 );
imageplot( convert_distance_color(D1,M,1), 'Update distance', 1,3,2 );
imageplot( convert_distance_color(min(D,D1),M,1), 'New distance', 1,3,3 );
% update
D = min(D,D1);


%EXO
%% Iterate the sampling process to add more and more points.
nmax = 200;
ndisp = [20 50 100 200];
vertex = [1;1];
[D,Z,Q] = perform_fast_marching(1./W, vertex);
clf; u = 1;
vertex_svg = {};
faces_svg = {};
for k=2:nmax
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
        u = u+1;  
        % compute the Delaunay triangulation
        [D1,Z,Q] = perform_fast_marching(1./W, vertex);
        vertex_svg{end+1} = vertex;
        faces_svg{end+1} = compute_voronoi_triangulation(Q,vertex);
    end
end
%EXO

%EXO
%% Display the geodesic Delaunay triangulation corresponding to the
%% sampling
clf;
for i=1:4
    subplot(2,2,i);
    plot_triangulation(vertex_svg{i},faces_svg{i}, M);
end
colormap jet(256);
%EXO

%% Lloyd Relaxation
% The farthest point sampling strategy is greedy and does not move the
% positions of the points once they are seeded.

%%
% To enhance the sampling, it is possible to relocate iteratively the
% points at the center of the Voronoi cells. This corresponds to the Lloyd
% algorithm, first developped for vector quantization.

%%
% First with a constant metric.

n = 512;
W = ones(n);

%%
% Seed random points.

p = 40;
vertex = floor(rand(2,p)*n-1)+1;

%%
% Compute Voronoi partition.

[D,Z,Q] = perform_fast_marching(1./W, vertex);


%%
% Display Vornoi cells.

clf; hold on;
imageplot(Q');
h = plot(vertex(1,:), vertex(2,:), 'k.');
set(h, 'MarkerSize', 15);
colormap(jet(256));

%%
% Re-center each point at the barycenter of its cell.

for i=1:p
    [x,y] = ind2sub(size(W), find(Q==i));
    vertex(:,i) = [mean(x);mean(y)];
end

%%
% Display updated partitions.

[D,Z,Q] = perform_fast_marching(1./W, vertex);
clf; hold on;
imageplot(Q');
h = plot(vertex(1,:), vertex(2,:), 'k.');
set(h, 'MarkerSize', 15);
colormap(jet(256));

%EXO
%% Perform the Lloyd iterative algorithm.
% metric
n = 512; W = ones(n);
% Seed random points.
p = 40;
vertex = floor(rand(2,p)*(n-1))+1;
disp_list = [1,2,5,20]; k = 1;
clf; hold on;
for i=1:max(disp_list)
    % Compute Voronoi partition.
    [D,Z,Q] = perform_fast_marching(1./W, vertex);
    if i==disp_list(k)
        subplot(2,2,k);
        imageplot(Q'); hold on;
        h = plot(vertex(1,:), vertex(2,:), 'k.');
        set(h, 'MarkerSize', 15);
        colormap(jet(256));
        k = k+1;
    end
    % Re-center each point at the barycenter of its cell.    
    for i=1:p
        [x,y] = ind2sub(size(W), find(Q==i));
        vertex(:,i) = [mean(x);mean(y)];
    end
end
%EXO

%CMT
% metric
n = 512; W = ones(n);
% Seed random points.
p = 150;
rand('state',1234);
vertex = floor(rand(2,p)*(n-1))+1;
disp_list = [1,2,5,30]; k = 1;
for i=1:max(disp_list)
    % Compute Voronoi partition.
    [D,Z,Q] = perform_fast_marching(1./W, vertex);
    if i==disp_list(k)
        clf;
        imageplot(Q'); hold on;
        h = plot(vertex(1,:), vertex(2,:), 'k.');
        set(h, 'MarkerSize', 15);
        colormap(jet(256));
        saveas(gcf, [rep 'lloyd-constant-' num2str(i) '.eps'], 'epsc');
        k = k+1;
    end
    % Re-center each point at the barycenter of its cell.    
    for i=1:p
        [x,y] = ind2sub(size(W), find(Q==i));
        vertex(:,i) = [mean(x);mean(y)];
    end
end
%CMT

%%
% Now we define a non-constent metric.

n = 256;
x = linspace(-1,1,n);
[Y,X] = meshgrid(x,x);
sigma = .3;
W = exp( -(X.^2+Y.^2)/(2*sigma^2) );
W = rescale(W,.02,1);

%%
% Display.

clf;
imageplot(W);

%EXO
%% Perform the Lloyd iterative algorithm.
% Seed random points.
p = 150;
vertex = floor(rand(2,p)*(n-1))+1;
disp_list = [1,2,5,30]; k = 1;
clf; hold on;
for i=1:max(disp_list)
    % Compute Voronoi partition.
    [D,Z,Q] = perform_fast_marching(1./W, vertex);
    if i==disp_list(k)
        subplot(2,2,k);
        imageplot(Q'); hold on;
        h = plot(vertex(1,:), vertex(2,:), 'k.');
        set(h, 'MarkerSize', 15);
        colormap(jet(256));
        k = k+1;
    end
    % Re-center each point at the barycenter of its cell.    
    for i=1:p
        [x,y] = ind2sub(size(W), find(Q==i));
        vertex(:,i) = [mean(x);mean(y)];
    end
end
%EXO


%CMT
% Seed random points.
p = 150;
rand('state',1234);
vertex = floor(rand(2,p)*(n-1))+1;
disp_list = [1,2,5,30]; k = 1;
for i=1:max(disp_list)
    % Compute Voronoi partition.
    [D,Z,Q] = perform_fast_marching(1./W, vertex);
    if i==disp_list(k)
        clf;
        imageplot(Q'); hold on;
        h = plot(vertex(1,:), vertex(2,:), 'k.');
        set(h, 'MarkerSize', 15);
        colormap(jet(256));
        saveas(gcf, [rep 'lloyd-gaussian-' num2str(i) '.eps'], 'epsc');
        k = k+1;
    end
    % Re-center each point at the barycenter of its cell.    
    for i=1:p
        [x,y] = ind2sub(size(W), find(Q==i));
        vertex(:,i) = [mean(x);mean(y)];
    end
end
%CMT