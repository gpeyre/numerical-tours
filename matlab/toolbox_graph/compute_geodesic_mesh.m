function [path,vlist,plist] = compute_geodesic_mesh(D, vertex, face, x, options)

% compute_geodesic_mesh - extract a discrete geodesic on a mesh
%
%   [path,vlist,plist] = compute_geodesic_mesh(D, vertex, face, x, options);
%
%   D is the set of geodesic distances.
%
%   path is a 3D curve that is the shortest path starting at x.
%   You can force to use a fully discrete descent using
%   options.method='discrete'.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
verb = getoptions(options, 'verb', 0);

if length(x)>1
    path = {}; vlist = {}; plist = {};
    for i=1:length(x)
        if length(x)>5     
            if verb
                progressbar(i,length(x));  
            end
        end
        [path{i},vlist{i},plist{i}] = compute_geodesic_mesh(D, vertex, face, x(i), options);
    end
    return;
end

method = getoptions(options, 'method', 'continuous');

[vertex,face] = check_face_vertex(vertex,face);

if strcmp(method, 'discrete')
    if isfield(options, 'v2v')
        vring = options.v2v;
    else
        vring =  compute_vertex_ring(face);
    end
    % path purely on edges
    vlist = x;
    vprev = D(x);
    while true
        x0 = vlist(end);
        r = vring{x0};
        [v,J] = min(D(r));
        x = r(J);
        if v>=vprev || v==0
            break;
        end
        vprev = v;
        vlist(end+1) = x;
    end
    path = vertex(:,vlist);
    plist = vlist*0+1;
    return;
end

%%% gradient descent on edges
% retrieve adjacency lists
m = size(face,2); n = size(vertex,2);
% precompute the adjacency datasets
if isfield(options, 'e2f')
    e2f = options.e2f;
else
    e2f = compute_edge_face_ring(face);
end
if isfield(options, 'v2v')
    v2v = options.v2v;
else
    v2v = compute_vertex_ring(face);
end
% initialize the paths
[w,f] = vertex_stepping(face, x, e2f, v2v, D);

vlist = [x;w]; 
plist = [1];
Dprev = D(x);

while true;
    % current triangle
    i = vlist(1,end);
    j = vlist(2,end);
    k = get_vertex_face(face(:,f),i,j);
    a = D(i); b = D(j); c = D(k);
    % adjacent faces
    f1 = get_face_face(e2f, f, i,k);
    f2 = get_face_face(e2f, f, j,k);
    % compute gradient in local coordinates
    x = plist(end); y = 1-x;
    gr = [a-c;b-c];
    % twist the gradient
    u = vertex(:,i) - vertex(:,k);
    v = vertex(:,j) - vertex(:,k);
    A = [u v]; A = (A'*A)^(-1);
    gr = A*gr;
    nx = gr(1); ny = gr(2);
    % compute intersection point
    etas = -y/ny;
    etat = -x/nx;
    s = x + etas*nx;
    t = y + etat*ny;
    if etas<0 && s>=0 && s<=1 && f1>0
        %%% CASE 1 %%%
        plist(end+1) = s;
        vlist(:,end+1) = [i k];
        % next face
        f = f1;
    elseif etat<0 && t>=0 && t<=1 && f2>0
        %%% CASE 2 %%%
        plist(end+1) = t;
        vlist(:,end+1) = [j k];
        % next face
        f = f2;
    else
        %%% CASE 3 %%%
        if a<=b
            z = i;            
        else
            z = j;
        end
        [w,f] = vertex_stepping( face, z, e2f, v2v, D);
        vlist(:,end+1) = [z w];  
        plist(end+1) = 1;   
    end
    Dnew = D(vlist(1,end))*plist(end) + D(vlist(2,end))*(1-plist(end));
    if Dnew==0 || (Dprev==Dnew && length(plist)>1)
        break;
    end
    Dprev=Dnew;
end

v1 = vertex(:,vlist(1,:));
v2 = vertex(:,vlist(2,:));

path = v1.*repmat(plist, [3 1]) + v2.*repmat(1-plist, [3 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w,f] = vertex_stepping(face, v, e2f, v2v, D)

% adjacent vertex with minimum distance
[tmp,I] = min( D(v2v{v}) ); w = v2v{v}(I);
f1 = e2f(v,w);
f2 = e2f(w,v);
if f1<0
    f = f2; return;
end
if f2<0
    f = f1; return;
end
z1 = get_vertex_face(face(:,f1),v,w);
z2 = get_vertex_face(face(:,f2),v,w);
if D(z1)<D(z2);
    f = f1; 
else
    f = f2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = get_vertex_face(f,v1,v2)

if nargin==2
    v2 = v1(2); v1 = v1(1);
end
k = setdiff(f, [v1 v2]);
if length(k)~=1
    error('Error in get_vertex_face');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = get_face_face(e2f, f, i,j)

f1 = e2f(i,j); f2 = e2f(j,i);
if f==f1
    f = f2;
else
    f = f1;
end



