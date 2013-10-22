function W = compute_mesh_weight(vertex,face,type,options)

% compute_mesh_weight - compute a weight matrix
%
%   W = compute_mesh_weight(vertex,face,type,options);
%
%   W is sparse weight matrix and W(i,j)=0 is vertex i and vertex j are not
%   connected in the mesh.
%
%   type is either 
%       'combinatorial': W(i,j)=1 is vertex i is conntected to vertex j.
%       'distance': W(i,j) = 1/d_ij^2 where d_ij is distance between vertex
%           i and j.
%       'conformal': W(i,j) = cot(alpha_ij)+cot(beta_ij) where alpha_ij and
%           beta_ij are the adjacent angle to edge (i,j)
%
%   If options.normalize=1, the the rows of W are normalize to sum to 1.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
[vertex,face] = check_face_vertex(vertex,face);

nface = size(face,1);
n = max(max(face));

if isfield(options, 'verb')
    verb = options.verb;
else
    verb = n>5000;
end

if nargin<3
    type = 'conformal';
end

switch lower(type)
    case 'combinatorial'
        W = triangulation2adjacency(face);
    case 'distance'
        W = my_euclidean_distance(triangulation2adjacency(face),vertex);
        W(W>0) = 1./W(W>0);
        W = (W+W')/2; 
    case 'conformal'
        % conformal laplacian
        W = sparse(n,n);
        ring = compute_vertex_face_ring(face);
        for i = 1:n
            if verb
                progressbar(i,n);
            end
            for b = ring{i}
                % b is a face adjacent to a
                bf = face(:,b);
                % compute complementary vertices
                if bf(1)==i
                    v = bf(2:3);
                elseif bf(2)==i
                    v = bf([1 3]);
                elseif bf(3)==i
                    v = bf(1:2);
                else
                    error('Problem in face ring.');
                end
                j = v(1); k = v(2);
                vi = vertex(:,i);
                vj = vertex(:,j);
                vk = vertex(:,k);
                % angles
                alpha = myangle(vk-vi,vk-vj);
                beta = myangle(vj-vi,vj-vk);
                % add weight
                W(i,j) = W(i,j) + cot( alpha );
                W(i,k) = W(i,k) + cot( beta );
            end
        end
    otherwise
        error('Unknown type.')
end

if isfield(options, 'normalize') && options.normalize==1
    W = diag(sum(W,2).^(-1)) * W;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta = myangle(u,v);

du = sqrt( sum(u.^2) );
dv = sqrt( sum(v.^2) );
du = max(du,eps); dv = max(dv,eps);
beta = acos( sum(u.*v) / (du*dv) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = my_euclidean_distance(A,vertex)

if size(vertex,1)<size(vertex,2)
    vertex = vertex';
end

[i,j,s] = find(sparse(A));
d = sum( (vertex(i,:) - vertex(j,:)).^2, 2);
W = sparse(i,j,d);  