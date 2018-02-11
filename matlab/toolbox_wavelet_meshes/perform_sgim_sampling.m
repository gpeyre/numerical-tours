function gim = perform_sgim_sampling(signal_original, pos_sphere, face, n, options)

% perform_sgim_sampling - generate a geometry image from a spherical parameterization
%
%   gim = perform_sgim_sampling(signal_original, pos_sphere, face, n, options);
%
%   'signal_original' are typically vertex location of the original mesh, but it can 
%   also be any information you want to re-sample on a regular grid (normal, etc).
%   'pos_sphere' are location on the sphere for these points.
%   'face' is the face data structure of the mesh.
%   'n' is width of the GIM.
%
%   Uses the spherical geometry image datastructure introduced in
%       Spherical parametrization and remeshing.
%       E. Praun, H. Hoppe.
%       SIGGRAPH 2003, 340-349.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<4
	n = 64;
end

options.null = 0;
sampling_type = getoptions(options, 'sampling_type', 'area');


if size(signal_original,2)~=size(pos_sphere,2)
    error('Original and spherical meshes must be of same size.');
end

% compute sampling location on the image
% disp('Computing planar sampling locations.');
posw = perform_spherial_planar_sampling(pos_sphere, sampling_type, options);

% perform 4-fold symmetry 
sym = { [0,1], [0,-1], [1, 0], [-1, 0] };
for i=1:length(sym)
    c = sym{i};
    if c(1)==0
        I = find(posw(2,:)*sign(c(2))>=0);
    else
        I = find(posw(1,:)*sign(c(1))>=0);
    end
    posi = posw(:,I);
    posi = perform_symmetry(posi,c);
    posw = [posw, posi];
    signal_original = [signal_original, signal_original(:,I)];
end

% crop a bit to speed up ..
m = 1.5;
I = find( posw(1,:)<m & posw(1,:)>-m & posw(2,:)<m & posw(2,:)>-m);
posw = posw(:,I);
signal_original = signal_original(:,I);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLATE THE CENTER (orignal triangulation)
% remove the faces that cross the boundary
edges = [face(1:2,:), face(2:3,:), face([1 3],:)];
pos = pos_sphere;
J = find( pos(3,edges(1,:))<=0 & pos(3,edges(2,:))<=0 &     ...
          ( pos(1,edges(1,:)).*pos(1,edges(2,:))<=0   |     ...
            pos(2,edges(1,:)).*pos(2,edges(2,:))<=0 )  );
% find the corresponding faces
J = unique( mod(J-1,size(face,2))+1 );
% find the complementary
x = zeros(size(face,2),1); x(J) = 1; I = find(x==0);
% retrive face number
face1 = face; face1 = face1(:,I);
% compute the interpolation on this sub-set using original triangulation
posn = (n-1)*(posw+1)/2+1;      % sampling location in [1,n]²
gim1 = zeros( n, n, size(signal_original,1) );
% fprintf('Griding using original triangulation ');
for i=1:size(signal_original, 1)
    % fprintf('.');
    gim1(:,:,i) = griddata_arbitrary( face1, posn, signal_original(i,:), n, options );
end


% remove doublons
if 1
[tmp,I] = unique(posw(1,:)+pi*posw(2,:));
posw = posw(:,I);
signal_original = signal_original(:,I);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLATE THE BOUNDARY (delaunay triangulation)
% sampling location in [-1,1]
x = -1:2/(n-1):1;
[Y,X] = meshgrid(x,x);
% interpolate location
gim2 = zeros( n, n, size(signal_original,1) );
% fprintf('\nGriding using Delaunay triangulation ');
for i=1:size(signal_original, 1)
    % fprintf('.');
    gim2(:,:,i) = griddata( posw(1,:), posw(2,:), signal_original(i,:), X, Y );
end
% fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIX THE TWO
gim = gim1;
I = find( isnan(gim) );
gim(I) = gim2(I);

I = find( isnan(gim) );
gim(I) = 0;

% to keep delaunay uncomment this
% gim = gim2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% translate/scale to fit in a box [0,1]^3

for i=1:size(signal_original, 1)
    x = gim(:,:,i); x = x(:);
    m = min( x );
    gim(:,:,i) = gim(:,:,i) - m;
end
gim = rescale(gim);


function y = perform_symmetry(x,c)
% y = 2*c - x
y(1,:) = 2*c(1) - x(1,:);
y(2,:) = 2*c(2) - x(2,:);

function M = griddata_arbitrary(face,vertex,v,n, options)

% griddata_arbitrary - perform interpolation of a triangulation on a regular grid
%
%   M = griddata_arbitrary(face,vertex,v,n);
%
%   n is the size of the image
%   vertex is a (2,n) list of points, they
%       are assumed to lie in [1,n]^2.
%   face is a (3,p) list of triangular faces.
%
%   If options.remove_nan==1, remove NaN non interpolated values using a
%   heat diffusion.
%   
%   Copyright (c) 2004 Gabriel Peyrée

options.null = 0;
verb = getoptions(options, 'verb', 1);
remove_nan = getoptions(options, 'remove_nan', 1);

if nargin<4
    error('4 arguments required.');
end
if size(face,2)==3 && size(face,1)~=3
    face = face';
end
if size(face,1)~=3
    error('Works only with triangulations.');
end
if size(vertex,2)==2 && size(vertex,1)~=2
    vertex = vertex';
end
if size(vertex,1)~=2
    error('Works only with 2D triangulations.');
end

if min(vertex(:))>=0 && max(vertex(:))<=1
    % date in [0,1], makes conversion
    vertex = vertex*(n-1)+1;
end

nface = size(face,2);
face_sampling = 0;
if length(v)==nface
    face_sampling = 1;
end

M = zeros(n);
Mnb = zeros(n);
for i=1:nface
    if verb
        progressbar(i,nface);
    end
    T = face(:,i);          % current triangles
    P = vertex(:,T);        % current points
    V = v(T);               % current values
    % range
    selx = min(floor(P(1,:))):max(ceil(P(1,:)));
    sely = min(floor(P(2,:))):max(ceil(P(2,:)));
    % grid locations
    [Y,X] = meshgrid(sely,selx);
    pos = [X(:)'; Y(:)'];
    p = size(pos,2);    % number of poitns
    % solve barycentric coords
    warning off;
    c = [1 1 1; P]\[ones(1,p);pos];
    warning on;
    % find points inside the triangle
    I = find( c(1,:)>=-10*eps & c(2,:)>=-10*eps & c(3,:)>=-10*eps );
    pos = pos(:,I);
    c = c(:,I);
    % restrict to point inside the image
    I = find( pos(1,:)<=n & pos(1,:)>0 & pos(2,:)<=n & pos(2,:)>0 );
    pos = pos(:,I);
    c = c(:,I);
    % convert sampling location to index in the array
    J = sub2ind(size(M), pos(1,:),pos(2,:) );
    % assign values
    if ~isempty(J)
        if face_sampling
            M(J) = v(i);
        else
            M(J) = M(J) + V(1)*c(1,:) + V(2)*c(2,:) + V(3)*c(3,:);
        end
        Mnb(J) = Mnb(J)+1;
    end
end
I = find(Mnb>0);
M(I) = M(I)./Mnb(I);
I = find(Mnb==0);
M(I) = NaN;

if remove_nan
    Mask = isnan(M);
    M(Mask==1) = mean(M(Mask==0));
    M0 = M;
    h = ones(3)/9;
    niter = 20;
    for i=1:niter
        M = perform_convolution(M,h);
        M(Mask==0) = M0(Mask==0);
    end
end
