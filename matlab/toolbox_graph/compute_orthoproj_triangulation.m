function v = compute_orthoproj_triangulation(vertex, face, M, options)

% compute_orthoproj_triangulation - compute best approximation on a triangulation
%
%	v = compute_orthoproj_triangulation(vertex, face, M, options);
%	
%	v are the values on the vertex of the triangulation that generates the 
%	best linear approximation of the image M using the triangulation.
%
%   To compute only an interpolation (and not a L2 projection), set
%   options.interpolate = 1.
%
%	Copyright (c) 2008 Gabriel Peyre

options.null = 0;
verb = getoptions(options, 'verb', 0);
interpolate = getoptions(options, 'interpolate', 0);

if interpolate
    v = interp2(M, vertex(2,:), vertex(1,:));
    return;
end

nface = size(face,2);
nvert = size(vertex,2);
n = size(M,1);

% mass matrix
U = spalloc(nvert,nvert,15*nvert);
% inner product matrix
y = zeros(nvert,1);
tag = zeros(nvert,1);
for i=1:nface
    
    if verb
        progressbar(i,nface);
    end

    T = face(:,i);          % current triangles
    P = vertex(:,T);        % current points

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
    if not(isempty(J))
        for j=1:3        
        j1 = mod(j,3)+1;
        % inner product with itself
        ip = sum( c(j,:).^2);
        U(T(j),T(j)) = U(T(j),T(j)) + ip;
        % inner product beteen two basis functions        
        ip = sum( c(j,:).*c(j1,:));
        U(T(j),T(j1)) = U(T(j),T(j1)) + ip;
        U(T(j1),T(j)) = U(T(j1),T(j)) + ip;        
        % image inner product
        ipM = sum( M(J).*c(j,:) );
        y(T(j)) = y(T(j)) + ipM;
        tag(T(j)) = 1;
        end
    end
end
% v = U(tag==1,tag==1)\y(tag==1);
warning off;
v = U \ y;
warning on;

% remove NaN values using interpolation
if sum(isnan(v))>0
    v(isnan(v)) = interp2(M, vertex(2,isnan(v)), vertex(1,isnan(v)));
end
