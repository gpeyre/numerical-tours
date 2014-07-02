function ic = check_incircle_edge(vertex, face, edge)

% check_incicle_edge - compute "empty circle" property for a set of edges
%
%   ic = check_incicle_edge(vertex,face, edge);
%
%   ic(i)==1 if edge(:,i) is delaunay valid (boundary or empty circles or non convex).
%   It thus should be flipped if ic(i)==0.
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin<3
    edge = compute_edges(face);
end

n = size(edge, 2);
e2f = compute_edge_face_ring(face);

m = size(e2f,1);
f1 = full( e2f(edge(1,:)+(edge(2,:)-1)*m) );
f2 = full( e2f(edge(2,:)+(edge(1,:)-1)*m) );

if not( isempty( find(f1==0 & f2==0) ) )
    warning('Problem with triangulation');
end

ic = ones(n,1);
isconvex = ones(n,1);

I = find(f1>0 & f2>0);

f1 = f1(I); f2 = f2(I);

f1 = face(:,f1);
f2 = face(:,f2);

edge = edge(:,I);
n = length(I);

points1 = find_other(f2, edge(1,:), edge(2,:));
ic1 = check_incircle(vertex, f1, points1 );

points2 = find_other(f1, edge(1,:), edge(2,:));
ic2 = check_incircle(vertex, f2, points2 );

ic(I) = ic1 & ic2;

% check for convexity
v = vertex(:,points2) - vertex(:,points1);
u1 = vertex(:,edge(1,:)) - vertex(:,points1);
u2 = vertex(:,edge(2,:)) - vertex(:,points1);
A = zeros(2,2,n); B = zeros(2,2,n);
A(:,1,:) = reshape(v, [2 1 n]);
A(:,2,:) = reshape(u1, [2 1 n]);
B(:,1,:) = reshape(v, [2 1 n]);
B(:,2,:) = reshape(u2, [2 1 n]);
isconvex(I) = ( det3(A) .* det3(B) ) < 0;
% non convex should not be flipped
ic(isconvex==0) = 1;


%%
function f = find_other(f,a,b)

u = f - repmat(a, [3 1]); f(u==0) = 0;
u = f - repmat(b, [3 1]); f(u==0) = 0;
if sum( sum(f==0)~=2 )>0
    error('Problem with triangulation');
end
f = sum(f);