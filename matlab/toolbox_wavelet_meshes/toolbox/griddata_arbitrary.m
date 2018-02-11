function M = griddata_arbitrary(face,vertex,v,n, options)

% griddata_arbitrary - perform interpolation of a triangulation on a regular grid
%
%   M = griddata_arbitrary(face,vertex,v,n);
%
% 'n' is the size of the image
% 'vertex' are assumed to lie in [1,n]^2
%
%   Copyright (c) 2004 Gabriel Peyré

options.null = 0;
verb = getoptions(options, 'verb', 1);

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
    c = [1 1 1; P]\[ones(1,p);pos];
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
