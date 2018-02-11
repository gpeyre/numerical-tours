%% Geodesic Segmentation
% This tour explores the use of Fast Marching methods for image segmentation.


perform_toolbox_installation('signal', 'general', 'graph');

%CMT
rep = 'results/fastmarching_2d/';
if not(exist(rep))
    mkdir(rep);
end
%CMT


%% Segmentation Using Geodesic Ball
% It is possible to extract an object by growing a geodesic ball.

%%
% First we load an image.

n = 256;
name = 'cortex';
M = rescale( sum(load_image(name,n),3) );

%%
% Display.

clf;
imageplot(M);

%% 
% Starting point of the grodesic ball.

pstart = [154;175];

%%
% Choose a metric that is minimal for value of the image close to pstart.

W = abs(M-M(pstart(1),pstart(2)));
W = rescale( max(W,0.03), 0.01,1).^2;

%CMT
imwrite(rescale(W), [rep name '-propaseg-metric.png'], 'png');
%CMT

%%
% Compute the Fast Marching from the center.

clear options;
options.nb_iter_max = Inf;
options.end_points = [];
[D,S,Q] = perform_fast_marching(1./W, pstart, options);

%EXO
%% Display geodesic balls {x \ M(x)<T} for various T.
v = sort(D(:));
Tlist = v( round([.05 .1 .15 .25]*n^2) );
clf;
for i=1:4
    T = Tlist(i);
    A = repmat(M,[1 1 3]);
    I = find(D<T);
    A(I) = 1;
    A([I+n^2;I+2*n^2]) = 0;
    subplot(2,2,i);
    hold on;
    imageplot(A);
    [c,h] = contour(D<T, [.5 .5], 'b');
    set(h, 'LineWidth', 2.5); 
    axis('ij');
end
%EXO

%CMT
for i=1:4
    T = Tlist(i);
    A = repmat(M,[1 1 3]);
    I = find(D<T);
    A(I) = 1;
    A([I+n^2;I+2*n^2]) = 0;
    clf;
    hold on;
    imageplot(A);
    [c,h] = contour(D<T, [.5 .5], 'b');
    set(h, 'LineWidth', 2.5); 
    axis('ij');
    saveas(gcf, [rep name '-propaseg-' num2str(i), '.eps'],'epsc');
end
%CMT

%% Segmentation Using Voronoi Diagrams
% It is possible to perform the segmentation by using an edge stopping
% metric, and Vornoi diagram for several seeds.

%%
% Magnitude of the gradient.

mu = 2;
d = sqrt( sum( grad(perform_blurring(M,mu)).^2, 3) ); 
d = perform_blurring(d,mu);

%% 
% Edge stopping metric.

W = rescale( min(d,0.15), 0.01,1).^2;

%%
% Display the metric.

clf;
imageplot(W);

%%
% Starting points.

pstart = [[30;30] [139;86] [158;170] [128;134] [124;122]];


%%
% Perform propagation.

options.nb_iter_max = Inf;
options.end_points = [];
[D,S,Q] = perform_fast_marching(1./W, pstart, options);

%%
% Display Voronoi diagrams.

clf;
imageplot(Q);
colormap(jet(256));


%EXO
%% Display the level sets.
col = {'r' 'g' 'b' 'c' 'y' 'k' 'r:'};
clf;
hold on;
imageplot(M);
Vlist = unique(Q(:))';
for i=Vlist
    U = zeros(n); U(Q==i)=1;
    [c,h] = contour(U, [.5 .5], col{i});
    set(h, 'LineWidth', 4); 
end
axis ij;
hold off;
%EXO

%CMT
saveas(gcf, [rep name '-voronoi-seg-image.eps'], 'epsc');
%CMT