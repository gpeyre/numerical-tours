function path = compute_geodesic(D, x, options)

% compute_geodesic - extract a discrete geodesic in 2D and 3D
%
%   p = compute_geodesic(D,x,options);
%
%   D is the distance map.
%   x is the starting point.
%   p is the shortest path between x and the starting point y (wich should
%   satisfy D(y)=0).
%
%   Set options.method = 'discrete' if you want to use a pure discrete
%   gradient descent.
%
%   Copyright (c) 2007 Gabriel Peyre


options.null = 0;
method = getoptions(options, 'method', 'continuous');

if strcmp(method, 'discrete')
    path = compute_discrete_geodesic(D,x);
    return;
end

if size(x,1)>1 && size(x,2)>1
    % several geodesics
    if size(x,1)>3
        x = x'; % probably user mistake of dimensions
    end
    path = {};
    for i=1:size(x,2)
        path{end+1} = compute_geodesic(D, x(:,i), options);
    end
    return;
end

if size(D,3)>1
    % 3D
    path = extract_path_3d(D,x,options);
else
    % 2D
    path = extract_path_2d(D,x,options);
end

if size(path,1)>size(path,2)
    path = path';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function path = compute_discrete_geodesic(D,x)

% compute_discrete_geodesic - extract a discrete geodesic in 2D and 3D
%
%   path = compute_discrete_geodesic(D,x);
%
%   Same as extract_path_xd but less precise and more robust.
%
%   Copyright (c) 2007 Gabriel Peyre

nd = 2;
if size(D,3)>1
    nd = 3;
end

x = x(:);
path = round(x(1:nd));

% admissible moves
if nd==2
    dx = [1 -1 0 0];
    dy = [0 0 1 -1];
    dz = [0 0 0 0];
    d = cat(1,dx,dy);
    vprev = D(x(1),x(2));
else
    dx = [1 -1 0 0 0 0];
    dy = [0 0 1 -1 0 0];
    dz = [0 0 0 0 1 -1];
    d = cat(1,dx,dy,dz);
    vprev = D(x(1),x(2),x(3));
end

s = size(D);
while true
    x0 = path(:,end);
    x = repmat(x0,1,size(d,2))+d;
    if nd==2
        I = find(x(1,:)>0 & x(2,:)>0 & x(1,:)<=s(1) & x(2,:)<=s(2) );    
    else
        I = find(x(1,:)>0 & x(2,:)>0 & x(3,:)>0 & x(1,:)<=s(1) & x(2,:)<=s(2) & x(3,:)<=s(3) );
    end
    x = x(:,I);
    if nd==2
        I = x(1,:) + (x(2,:)-1)*s(1);
    else
        I = x(1,:) + (x(2,:)-1)*s(1) + (x(3,:)-1)*s(1)*s(2);
    end
    [v,J] = min(D(I));
    x = x(:,J);
    if v>vprev
        return;
    end
    vprev = v;
    path(:,end+1) = x;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function path = extract_path_2d(A,end_points,options)

% extract_path_2d - extract the shortest path using 
%   a gradient descent.
%
%   path = extract_path_2d(D,end_point,options);
%
%   'D' is the distance function.
%   'end_point' is ending point (should be integer). 
%
%   Copyright (c) 2004 Gabriel Peyré

options.null = 0;

if isfield(options, 'trim_path')
    trim_path = options.trim_path;
else
    trim_path = 1;
end

if isfield(options, 'stepsize')
    stepsize = options.stepsize;
else
    stepsize = 0.1;
end
if isfield(options, 'maxverts')
    maxverts = options.maxverts;
else
    maxverts = 10000;
end
str_options = [stepsize maxverts];

if size(end_points,1)~=2
    end_points = end_points';
end
if size(end_points,1)~=2
    error('end_points should be of size 2xk');
end

% gradient computation
I = find(A==Inf);
J = find(A~=Inf);
A1 = A; A1(I) = mmax(A(J));
global grad;
grad = compute_grad(A1);
grad = -perform_vf_normalization(grad);

% path extraction
path = stream2(grad(:,:,2),grad(:,:,1),end_points(2,:),end_points(1,:), str_options);
for i=1:length(path)
    path{i} = path{i}(:,2:-1:1);
end
if length(path)==1
    path = path{1};
end

if isfield(options, 'start_points')
    start_points = options.start_points;
else
    start_points = path(end,:);
end
start_points = start_points(:);

if 0 % trim_path
    % removing too verbose points
    d = compute_distance_to_points(path', start_points);
    % perform thresholding
    T = mmax(d)/300^2;
    I = find(d<T);
    if not(isempty(I))
        path = path(1:I(1), :);
        path = [path; start_points'];
    else
        path = path';
    end
end

% complete with a discrete extraction (nasty hack)
if size(path, 2)~=2 && size(path, 1)==2
    path = path';
end
path = [path; compute_discrete_geodesic(A, round(path(end,:)))'];

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        OLD METHOD USING ODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isfield(options, 'Tmax')
    Tmax = options.Tmax;
else
    Tmax = sum(size(A));
end

global start_points;
if isfield(options, 'start_points')
    start_points = options.start_points;
else
    start_points = [];
end

x = x(:);

% path is empty
path = x;


% n * mmax(A(J));
% x = xv;

options = odeset('Events',@event_callback);
[T,path] = ode113( @get_gradient, [0,Tmax], x, options);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction
function [value,isterminal,direction] = event_callback(t,y)

global start_points;
% compute distance to start points
d = compute_distance_to_points(y,start_points);
value = min(d);
if value<0.1
    value = -1;
end
isterminal = 1;
direction = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction
function g = get_gradient( t, y )

global grad;
gx = grad(:,:,1);
gy = grad(:,:,2);

n = length(gx);

% down/left corner of current cell
p = floor(y(1));
q = floor(y(2));
    
p = clamp(p,1,n-1);
q = clamp(q,1,n-1);

% residual in [0,1]
xx = y(1)-p;
yy = y(2)-q;

xx = clamp(xx,0,1);
yy = clamp(yy,0,1);

% compute gradient    
a = [gx(p,q), gy(p,q)];
b = [gx(p+1,q), gy(p+1,q)];
c = [gx(p,q+1), gy(p,q+1)];
d = [gx(p+1,q+1), gy(p+1,q+1)];     
g = ( a*(1-xx)+b*xx )*(1-yy) + ( c*(1-xx)+d*xx )*yy;
g = -g';






function path = extract_path_3d(A,end_points,options)

% extract_path_3d - extract the shortest path using 
%   a gradient descent.
%
%   path = extract_path_3d(A,x,options);
%
%   'A' is the distance function.
%   'x' is starting point (should be integer). 
%
%   Copyright (c) 2004 Gabriel Peyré


options.null = 0;

if isfield(options, 'trim_path')
    trim_path = options.trim_path;
else
    trim_path = 1;
end

% gradient computation
I = find(A==Inf);
J = find(A~=Inf);
A1 = A; A1(I) = mmax(A(J));
global gx;
global gy;
global gz;
[gy,gx,gz] = gradient(A1);

% normalize the gradient field
d = sqrt( gx.^2 + gy.^2 + gz.^2 );
I = find(d<eps);
d(I) = 1;
gx = gx./d; gy = gy./d; gz = gz./d;

% path extraction
options = [0.2 20000];
path = stream3(-gy,-gx,-gz,end_points(2,:),end_points(1,:),end_points(3,:), options);
for i=1:length(path)
     path{i} = path{i}(:,[2:-1:1 3]);
end
if length(path)==1
    path = path{1};
end

% test if path is long enough
v = path(end,:);
v1 = max(round(v),ones(1,3));
if( 0 && A1(v1(1),v1(2),v1(3))>0.1 )
    path1 = stream3(-gy,-gx,-gz,v(2),v(1),v(3), options);
    for i=1:length(path1)
        path1{i} = path1{i}(:,[2:-1:1 3]);
    end
    if length(path)>=1
        path1 = path1{1};
    end
    path = [path; path1];
end

if isfield(options, 'start_points')
    start_points = options.start_points;
else
    start_points = path(end,:);
end
start_points = start_points(:);

if trim_path
    % removing too verbose points
    d = compute_distance_to_points(path', start_points);
    % perform thresholding
    T = mmax(d)/300^2;
    I = find(d<T);
    if ~isempty(I)
        path = path(1:I(1), :);
        path = [path; start_points'];
    end
end

function v = mmax(x)
v = max(x(:));


function grad = compute_grad(M,options)

% compute_grad - compute the gradient of an image using central differences
%
% grad = compute_grad(M,options);
%
%   'options' is a structure:
%   - options.h is the sampling step size on both direction (default 1).
%   - options.h1 is the sampling step size on X direction (default 1).
%   - options.h2 is the sampling step size on Y direction (default 1).
%   - options.type is the kind of finite difference.
%       type==2 is fwd differences, ie.
%           y(i) = (x(i)-x(i-1))/h, with special
%           care at boundaries.
%       type==1 is forward differences bilinearly interpolated in the
%           middle of each pixel (be aware that you have a shift of 1/2 on X and Y for
%           the location of the gradient).
%       type==1 is backward differences bilinearly interpolated in the
%           middle of each pixel (be aware that you have a shift of -1/2 on X and Y for
%           the location of the gradient).
%
%   Copyright (c) 2004 Gabriel PeyrÈ

if nargin<2
    options.null = 0;
end

if ~isfield(options, 'h1')
    options.h1 = 1;
end
h1 = options.h1;
if ~isfield(options, 'h2')
    options.h2 = 1;
end
h2 = options.h2;

if isfield(options, 'h')
    h1 = options.h;
    h2 = options.h;
end

[n,p] = size(M);

if isfield(options, 'type')
    type = options.type;
else
    type = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% new code, use faster 2D differences
if type==1
    % central differences on X
    D1 = [M(2:end,:);M(end,:)];
    D2 = [M(1,:);M(1:end-1,:)];
    grad(:,:,1) = (D1-D2)/(2*h1); 
    grad(1,:,1) = ( 4*M(2,:) - 3*M(1,:) - M(3,:) )/(2*h1);
    grad(end,:,1) = -( 4*M(end-1,:) - 3*M(end,:) - M(end-2,:) )/(2*h1);
    % central differences on Y
    D1 = [M(:,2:end),M(:,end)];
    D2 = [M(:,1),M(:,1:end-1)];
    grad(:,:,2) = (D1-D2)/(2*h2); 
    grad(:,1,2) = ( 4*M(:,2) - 3*M(:,1) - M(:,3) )/(2*h2);
    grad(:,end,2) = -( 4*M(:,end-1) - 3*M(:,end) - M(:,end-2) )/(2*h2);
elseif type==2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % accumulate on Y
    MM = ( M + [M(:,2:end),M(:,end)] )/2;
    % fwd differences on X
    D1 = [MM(2:end,:);MM(end,:)];
    D2 = MM;
    grad(:,:,1) = (D1-D2)/h1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % accumulate on X
    MM = ( M + [M(2:end,:);M(end,:)] )/2;
    % fwd differences on Y
    D1 = [MM(:,2:end),MM(:,end)];
    D2 = MM;
    grad(:,:,2) = (D1-D2)/h2;
elseif type==3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % accumulate on Y
    MM = ( M + [M(:,1),M(:,1:end-1)] )/2;
    % fwd differences on X
    D1 = MM;
    D2 = [MM(1,:);MM(1:end-1,:)];
    grad(:,:,1) = (D1-D2)/h;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % accumulate on Y
    MM = ( M + [M(1,:);M(1:end-1,:)] )/2;
    % fwd differences on Y
    D1 = MM;
    D2 = [MM(:,1),MM(:,1:end-1)];
    grad(:,:,2) = (D1-D2)/h;
else
    error('This kind of differences is not supported.');
end
 
return;

if type~=1
    % compute the difference in the center of each square
    h = zeros(n,p,2);
    
    for j=1:p-1
        h(:,j,1) = ( grad(:,j,1)+grad(:,j+1,1) )/2;
    end
    for i=1:n-1
        h(i,:,2) = ( grad(i,:,2)+grad(i+1,:,2) )/2;
    end
    
    if type==2          % fwd differences
        grad(1:n-1,1:p-1,:) = h(1:n-1,1:p-1,:);        
    elseif type==3      % bwd differences
        grad(2:n,2:p,:) = h(2:n,2:p,:);        
    end   
    
end


function v2 = perform_vf_normalization(v1)

% perform_vf_normalization - renormalize a vf.
%
%   v2 = perform_vf_normalization(v1);
%
%   Copyright (c) 2004 Gabriel PeyrÈ


n = sum(v1.^2,3);
I = find(n<eps);
n(I) = 1;
v2 = v1./repmat( sqrt(n), [1 1 2] );