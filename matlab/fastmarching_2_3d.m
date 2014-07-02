%% Fast Marching in 3D
% This tour explores the use of Fast Marching methods in 2D.

perform_toolbox_installation('signal', 'general', 'graph');

%% 3D Volumetric Datasets
% A volumetric data is simply a 3D array.

%%
% We load a volumetric data.

name = 'vessels';
options.nbdims = 3;
M = read_bin(name, options);
M = rescale(M);

%% 
% Size of the image (here it is a cube).

n = size(M,1);

%%
% Such a volumetric dataset is more difficult to visualize than a standard 2D image. 
% You can render slices along each X/Y/Z direction.

clf;
imageplot(M(:,:,50), 'X/Y slice', 1, 3, 1);
imageplot(squeeze(M(:,50,:)), 'X/Z slice', 1, 3, 2);
imageplot(squeeze(M(50,:,:)), 'Y/Z slice', 1, 3, 3);

%%
% We can display some horizontal slices.

slices = round(linspace(10,n-10,4));
clf;
for i=1:length(slices)
    s = slices(i);
    imageplot( M(:,:,s), strcat(['Z=' num2str(s)]), 2,2,i );
end


%%
% You can also perform a volumetric rendering.
% In order to do so, you need to set up a correct alpha mapping to make transparent some parts of the volume. 
% Here, each time the options.center value is increased.

clf;
h = vol3d('cdata',M,'texture','2D');
view(3); axis off;
% set up a colormap
colormap bone(256);
% set up an alpha map
options.sigma = .08; % control the width of the non-transparent region
options.center = .4; % here a value in [0,1]
a = compute_alpha_map('gaussian', options); % you can plot(a) to see the alphamap
% refresh the rendering
vol3d(h);

%EXO 
%% Try with other alphamapping and colormapping
clf;
h = vol3d('cdata',M,'texture','2D');
view(3); axis off;
% set up a colormap
colormap bone(256);
% set up an alpha map
options.sigma = .08; % control the width of the non-transparent region
options.center = .6; % here a value in [0,1]
a = compute_alpha_map('gaussian', options); % you can plot(a) to see the alphamap
% refresh the rendering
vol3d(h);
%EXO 


%%
% We can display an isosurface of the dataset (here we sub-sample to speed
% up the computation).

sel = 1:2:n;
clf;
isosurface( M(sel,sel,sel), .5);
axis('off');


%% 3D Shortest Paths
% The definition of shortest path extend to any dimension, for instance 3D.
   
%%
% Geodesic distances can be computed on a 3D volume using the 
% Fast Marching. The important point here is to define the correct
% potential field W that should be large in the region where you want the front to move fast. 
% It means that geodesic will follow these regions. 


%%
% Select the starting points.

delta = 5;
start_point = [107;15;delta];

%%
% Compute the (inverse of the) potential that is small  close
% to the value of |M| at the selected point

W = abs( M - M(start_point(1),start_point(2),start_point(3)) );

%%
% Rescale above some threshold to avoid too small potentials.

W = rescale(W,1e-2,1);

%%
% Perform the front propagation.

options.nb_iter_max = Inf;
[D,S] = perform_fast_marching(1./W, start_point, options);

%%
% In order to extract a geodesic, we need to select an ending point
% and perform a descent of the distance function D from this starting point.
% The selection is done by choosing a point of low distance value in the slice D(:,:,n-delta).

clf;
imageplot(D(:,:,delta), '', 1,2,1);
imageplot(D(:,:,n-delta), '', 1,2,2);
colormap(jet(256));

%EXO
%% select the point (x,y) of minimum value in the slice |D(:,:,n-delta)|.
%% hint: use functions 'min' and 'ind2sub'
d = D(:,:,n-delta);
[tmp,I] = min(d(:));
[x,y] = ind2sub([n n],I(1));
end_point = [x;y;n-delta];
%EXO

%%
% Extract the geodesic by discrete descent.

options.method = 'discrete';
minpath = compute_geodesic(D,end_point,options);

%%
% Draw the path.

Dend = D(end_point(1),end_point(2),end_point(3));
D1 = double( D<=Dend );
clf;
plot_fast_marching_3d(M,D1,minpath,start_point,end_point);

%EXO
%% Select other starting points. In order to do so, ask the user to 
%% click on a starting point in a given horizontal slice |W(:,:,delta)|.
%% You can do this by using |ginput(1)|
%% on the plane |Z=delta|.
if 0
    clf; imageplot(M(:,:,delta));
    title('Pick starting point');
    start_point = round( ginput(1) );
    start_point = [start_point(2); start_point(1); delta];
end
%EXO