function faces = compute_voronoi_triangulation(Q, vertex)

% compute_voronoi_triangulation - compute a triangulation
%
%   face = compute_voronoi_triangulation(Q);
%
%   Q is a Voronoi partition function, computed using
%       perform_fast_marching.
%   face(:,i) is the ith face.
%
%   Works in 2D and in 3D.
%
%   Copyright (c) 2008 Gabriel Peyre

d = nb_dims(Q);

if d==2
    faces = compute_voronoi_triangulation_2d(Q, vertex);
elseif d==3
    faces = compute_voronoi_triangulation_3d(Q, vertex);
else
    error('Works only for 2D or 3D data.');
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function faces = compute_voronoi_triangulation_3d(Q, vertex)

[dx dy dz] = meshgrid(0:1,0:1,0:1);
V = [];
for i=1:8
    v = Q(1+dx(i):end-1+dx(i),1+dy(i):end-1+dy(i),1+dz(i):end-1+dz(i)); 
    V = [V v(:)];
end

V = sort(V,2);
V = unique(V, 'rows');
V = V( prod(V,2)>0 ,:);

d = V(:,1)*0;
for i=1:7
   d = d+(V(:,i)~=V(:,i+1));
end
if sum(d==4)>1
%    warning('Problem with triangulation.');
end

% split 5 folds into 2
I = find(d==4);
faces = [];
for i=1:length(I)
    v = unique(V(I(i),:));
    x = vertex(:,v); % points
    % barycenter
    a = sum( x, 2 ) / 5;
    x = x-repmat(a, [1 size(x,2)]);
    t = atan2(x(2,:),x(1,:));
    [tmp,s] = sort(t);
    faces = [faces v( s([1 2 3 5]))'];
    faces = [faces v( s([3 4 1 5]))'];
    
end
%    faces = [V(I,1:3);V(I,[3 4 1])];

% remaining triangles
V = V( d==3 ,:);
for i=1:size(V,1)
    V(i,1:4) = unique(V(i,:));
end
faces = [faces, V(:,1:4)'];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function faces = compute_voronoi_triangulation_2d(Q, vertex)

V = [];
v = Q(1:end-1,1:end-1); V = [V v(:)];
v = Q(2:end,1:end-1); V = [V v(:)];
v = Q(1:end-1,2:end); V = [V v(:)];
v = Q(2:end,2:end); V = [V v(:)];

V = sort(V,2);
V = unique(V, 'rows');
V = V( prod(V,2)>0 ,:);


d = (V(:,1)~=V(:,2)) + (V(:,2)~=V(:,3)) + (V(:,3)~=V(:,4));
if sum(d==3)>1
%    warning('Problem with triangulation.');
end

% split squares into 2 triangles
I = find(d==3);
faces = [];
for i=1:length(I)
    v = V(I(i),:);
    x = vertex(:,v); % points
    % barycenter
    a = sum( x, 2 ) / 4;
    x = x-repmat(a, [1 size(x,2)]);
    t = atan2(x(2,:),x(1,:));
    [tmp,s] = sort(t);
    faces = [faces; v( s([1 2 3]))];
    faces = [faces; v( s([3 4 1]))];
    
end
%    faces = [V(I,1:3);V(I,[3 4 1])];

% remaining triangles
V = V( d==2 ,:);
for i=1:size(V,1)
    V(i,1:3) = unique(V(i,:));
end
faces = [faces; V(:,1:3)];

% return in correct format
faces = faces';
