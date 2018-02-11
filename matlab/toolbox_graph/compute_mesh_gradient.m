function G = compute_mesh_gradient(vertex,face,type,options)

% compute_mesh_laplacian - compute a gradient matrix
%
%   G = compute_mesh_gradient(vertex,face,type,options);
%
%   G is an (m,n) matrix where n is the number of vertex and m the number
%   of edges in the mesh (we assume edges are oriented).
%
%   One has G((i,j),k)=sqrt(W(i,j)) if i==k and 
%   G((i,j),k)=-sqrt(W(i,j)) if j==k.
%   (in this definition we assume that i<j)
%
%   One has L=G'*G where L is the laplacian matrix.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;

if isfield(options, 'normalize')
    normalize = options.normalize;
else
    normalize = 1;
end

options.normalize = 0;
W = compute_mesh_weight(vertex,face,type,options);

%% compute list of edges
[i,j,s] = find(sparse(W));
I = find(i<j);
i = i(I);
j = j(I);
s = sqrt(s(I));
% number of edges
m = length(i);
% number of vertices
n = size(W,1);

%% build sparse matrix
s = [s; -s];
is = [(1:m)'; (1:m)'];
js = [i(:); j(:)];
G = sparse(is,js,s,m,n);

if normalize
    G = G*diag(sum(W,2).^(-1/2));
end