function [Q,DQ, voronoi_edges, edges_id, lambda] = compute_voronoi_mesh(vertex,faces, start_points, options)

% compute_voronoi_mesh - compute the exact voronoi segmentation of a 3D mesh
%
%   [Q,DQ, voronoi_edges, edges_id, edges_lambda] = compute_voronoi_mesh(vertex,faces,start_points, options);
%
%   Q(i) is the voronoi cell number of the i th point in the mesh.
%   DQ(:,1) is the distance to the closest point.
%   DQ(:,2) is the distance to the 2nd closest point.
%   DQ(:,3) is the distance to the 3rd closest point.
%
%   voronoi_edges is the of segment that define voronoi cells boundary.
%
%   See also perform_fast_marching_mesh.
%
%   Copyright (c) 2007 Gabriel Peyre

[vertex,faces] = check_face_vertex(vertex,faces);

n = size(vertex,2);
nstart = length(start_points);
options.null = 0;
verb = getoptions(options, 'verb', 0);

options.end_points = [];
[D0,S,Q0] = perform_fast_marching_mesh(vertex, faces, start_points, options);

D = [];
for i=1:nstart
    if verb
        progressbar(i,nstart);
    end
    options.dmax = max( D0(Q0==start_points(i)) )*1.1;
    [D(:,i),S,Q] = perform_fast_marching_mesh(vertex, faces, start_points(i), options);
end

% first hit
D1 = D;
Q = []; DQ = [];
for i=1:3
    [DQ(:,i),Q(:,i)] = min(D1,[],2);
    Qind = (1:n)' + (Q(:,i)-1)*n;
    D1(Qind) = Inf;
end
Q1 = Q(:,1); % first voronoi

edge = compute_edges(faces);
i = edge(1,:); j = edge(2,:);
% crossings
I = find(Q1(i)~=Q1(j));
i = i(I); j = j(I);
edges_id = [i(:)'; j(:)'];
qi = Q(i,1);  qj = Q(j,1); 
a = D( i' + (qi-1)*n );
b = D( j' + (qi-1)*n );
c = D( i' + (qj-1)*n );
d = D( j' + (qj-1)*n );
lambda = (c-a) ./ (b-a+c-d);
lambda = clamp(lambda(:)',0,1);
% lambda = 1/2+lambda*0;
% point locations
vcross = vertex(:,i) .* repmat(1-lambda, [3 1]) + vertex(:,j) .* repmat(lambda, [3 1]);
% convert edges to crossings number 
E2C = sparse(i,j,1:length(i));

A = (Q1(faces(1,:))~=Q1(faces(2,:))) + (Q1(faces(2,:))~=Q1(faces(3,:))) + (Q1(faces(1,:))~=Q1(faces(3,:)));
vcross1 = [];
vcross2 = [];
% find faces where crossings occurs
I = find(A==2);
for i=I(:)'
    e = [];
    for j=1:3
        j1 = mod(j,3)+1;
        if Q1(faces(j,i))~=Q1(faces(j1,i))
            q1 = min( faces(j,i), faces(j1,i) );
            q2 = max( faces(j,i), faces(j1,i) );
            ne = E2C( q1,q2 );
            e(:,end+1) = vcross(:,ne);
        end
    end
    vcross1(:,end+1) = e(:,1);
    vcross2(:,end+1) = e(:,2);
end
% find faces with triple crossings
I = find(A==3);
for i=I(:)'
    e = [];
    for j=1:3
        j1 = mod(j,3)+1;
        q1 = min( faces(j,i), faces(j1,i) );
        q2 = max( faces(j,i), faces(j1,i) );
        ne = E2C( q1,q2 );
        e(:,end+1) = vcross(:,ne);
    end
    ee = mean(e,2);
    vcross1(:,end+1:end+3) = e;
    vcross2(:,end+1:end+3) = [ee ee ee];
end

voronoi_edges = [vcross1; vcross2];

