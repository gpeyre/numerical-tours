%% Shape Retrieval with Geodesic Descriptors
% This numerical tour explores the use of geodesic distances within shapes
% to perform shape retrieval.

%%
% This tour is mostly inspired from the following work:

%%
% _Matching 2D and 3D Articulated Shapes using Eccentricity_, 
% A. Ion, N. M. Artner, G. Peyre, W. G. Kropatsch and L. Cohen,
% Preprint Hal-00365019, January 2009.

perform_toolbox_installation('signal', 'general', 'graph');

%CMT
rep = 'results/shapes_geodesic_descriptors/';
if not(exist(rep))
    mkdir(rep);
end
%CMT

%% Geodesic Distances Within a Binary Shape
% By restricting shortest path to lie within a shape, one create a geodesic
% metric that is different from the Euclidean one if the shape is not
% convex.

%%
% A binary shape is represented as a binary image.

n = 200;
name = 'centaur1';
M = load_image(name,n);
M = perform_blurring(M,5);
M = double( rescale( M )>.5 );
if M(1)==1 
    M = 1-M;
end


%%
% Display the shape.

clf;
imageplot(-M);

%CMT
imwrite(rescale(-M), [rep 'geodescr-' name '-original.png'], 'png');
%CMT

%%
% Compute its boundary

bound = compute_shape_boundary(M);
nbound = size(bound,2);

%% 
% Parameters for the Fast Marching: constant speed |W|, but retricted using |L| to the
% inside of the shape.

W = ones(n);
L = zeros(n)-Inf; L(M==1) = +Inf;


%%
% Initial point for the geodesic computation.

start_points = [95; 20];

%%
% Compute the geodesic distance without constraint using Fast Marching.
% It is simply the Euclidean distance.

options.constraint_map = [];
D0 = perform_fast_marching(W, start_points, options);
D0(M==0) = Inf;

%%
% Display Euclidean distance.

clf;
options.display_levelsets = 1;
options.pstart = start_points;
options.nbr_levelsets = 30;
display_shape_function(D0, options);

%CMT
saveas(gcf, [rep 'geodescr-' name '-disteucl.eps'], 'epsc');
%CMT

%%
% Compute the geodesic distance with constraints using Fast Marching.

options.constraint_map = L;
D = perform_fast_marching(W, start_points, options);

%%
% Display geodesic distance.

clf;
options.nbr_levelsets = 60;
display_shape_function(D, options);

%CMT
saveas(gcf, [rep 'geodescr-' name '-distgeod.eps'], 'epsc');
%CMT

%EXO
%% Using |options.nb_iter_max| display the progression of the Fast
%% Marching.
p = sum(M(:)==1);
qlist = ceil([.1 .3 .7 1]*p);
clf;
for i=1:4
    options.nb_iter_max = qlist(i);
    D = perform_fast_marching(W, start_points, options);
    subplot(2,2,i);
    display_shape_function(perform_hist_eq(D, 'linear')); 
    hold on; 
    h = plot(bound(2,:), bound(1,:), 'k');
    set(h, 'LineWidth', 2);
    axis('ij');
end
options.nb_iter_max = Inf;
%EXO

%CMT
p = sum(M(:)==1);
qlist = ceil([.1 .3 .7 1]*p);
clf;
for i=1:4
    options.nb_iter_max = qlist(i);
    D = perform_fast_marching(W, start_points, options);
    clf;
    display_shape_function(perform_hist_eq(D, 'linear')); 
    hold on; 
    h = plot(bound(2,:), bound(1,:), 'k');
    set(h, 'LineWidth', 2);
    axis('ij');
    saveas(gcf, [rep 'geodescr-' name '-fm-' num2str(i) '.eps'], 'epsc');
end
options.nb_iter_max = Inf;
%CMT

%%
% Compute a geodesic curve.

end_points = [27;112];
p = compute_geodesic(D,end_points);

%%
% Display the path.

ms = 30; lw = 3;
clf; hold on;
imageplot(1-M);
h = plot(end_points(2),end_points(1), '.b'); set(h, 'MarkerSize', ms);
h = plot(start_points(2),start_points(1), '.r'); set(h, 'MarkerSize', ms);
h = plot( p(2,:), p(1,:), 'g' ); set(h, 'LineWidth', lw);
axis ij;

%EXO
%% Compute curves joining the start point to several points along the
%% boundary.
npaths = 30;
sel = round(linspace(1,nbound+1,npaths+1)); sel(end) = [];
end_points = bound(:,sel);
clf; hold on;
imageplot(1-M);
for i=1:npaths
    p = compute_geodesic(D,end_points(:,i));
    h = plot( p(2,:), p(1,:), 'g' ); set(h, 'LineWidth', lw);
end
h = plot(start_points(2),start_points(1), '.r'); set(h, 'MarkerSize', ms);
h = plot(end_points(2,:),end_points(1,:), '.b'); set(h, 'MarkerSize', ms);
axis ij;
%EXO

%% Local Geodesic Descriptors
% In order to build shape signatures, we compute geodesic distance to all
% the points on the boundary. We then retrieve some caracteristics from
% these geodesic distance map.

%%
% Select a uniform set of points on the boundary.

nb_samples = 600;
sel = round(linspace(1,nbound+1,nb_samples+1)); sel(end) = [];
samples = bound(:,sel);


%EXO
%% Build a collection |E| of distance maps, so that |E(:,:,i)| is the
%% geodesic distance to |samples(:,i)|. 
E = zeros(n,n,nb_samples);
for i=1:nb_samples
    % progressbar(i,nb_samples);
    [d,tmp] = perform_fast_marching(W, samples(:,i), options);
    d(M==0) = 0;
    d(d==Inf) = 0;
    d(d>1e5) = 0;
    E(:,:,i) = d;
end
%EXO

%%
% normalize distances.

E = E/mean(E(:));


%%
% Display some locations

points = [[80;20] [95;112] [156;42]];
col = {'r', 'g', 'b', 'k'};
clf; hold on;
imageplot(-M);
for i=1:3
    h = plot(points(2,i), points(1,i), [col{i} '.']);
    set(h, 'MarkerSize', 40);
end
axis('ij');


%CMT
saveas(gcf, [rep 'geodescr-' name '-samples-loc.eps'], 'epsc');
%CMT

%%
% Display three different features at some locations.

clf;
col = {'r', 'g', 'b', 'k'};
for i=1:3
    subplot(3,1,i);    
    d = E(points(1,i),points(2,i), :); 
    u = hist(d(:), 15); axis tight;
    h = bar(u, col{i}); axis('tight');
    set(gca, 'XTickLabel', []);
end

%CMT
saveas(gcf, [rep 'geodescr-' name '-samples-descr.eps'], 'epsc');
%CMT


%% Global Geodesic Descriptors
% One can retain a single statistic from the local descriptors, such as the
% min, max, mean or median values. The histogram of these values are the
% global descriptors.

%%
% Compute several statistics.

clear A;
A{1} = max(E,[],3);
A{2} = min(E,[],3);
A{3} = mean(E,3);
A{4} = median(E,3);
titles = {'Max', 'Min', 'Mean', 'Median'};


%%
% Display as images.

nbr = [20 5 30 30];
options.pstart = [];
clf;
for i=1:4
    subplot(2,2,i);
    options.nbr_levelsets = nbr(i);
    display_shape_function(A{i}, options);
    title(titles{i});
end
colormap jet(256);

%CMT
clf;
for i=1:4
    options.nbr_levelsets = nbr(i);
    display_shape_function(A{i}, options);
    saveas(gcf, [rep 'geodescr-' name '-' titles{i} '.eps'], 'epsc');
end
%CMT

%%
% Display histograms of the statistics.

clf;
for i=1:4
    u = A{i}(M==1); u = u(u>0);
    subplot(4,1,i);
    hist(u, 40); axis('tight');
    title(titles{i});
end

%% Shape Retrieval using Geodesic Historams.
% One can use the histograms of Eccentricity for shape retrieval.

%EXO
%% Load a library of shapes. Compute the different histograms for these
%% shapes.
%EXO

%EXO
%% Perform the retrieval by comparing the histogram. Test diffetent metrics
%% for the retrieval.
%EXO
