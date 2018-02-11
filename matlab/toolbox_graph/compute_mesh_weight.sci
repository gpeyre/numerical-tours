function W = compute_mesh_weight(vertex,face,laptype,options)

// compute_mesh_weight - compute a weight matrix
//
//   W = compute_mesh_weight(vertex,face,laptype,options);
//
//   W is sparse weight matrix and W(i,j)=0 is vertex i and vertex j are not
//   connected in the mesh.
//
//   laptype is either 
//       'combinatorial': W(i,j)=1 is vertex i is conntected to vertex j.
//       'distance': W(i,j) = 1/d_ij^2 where d_ij is distance between vertex
//           i and j.
//       'conformal': W(i,j) = cot(alpha_ij)+cot(beta_ij) where alpha_ij and
//           beta_ij are the adjacent angle to edge (i,j)
//
//   If options.normalize=1, the the rows of W are normalize to sum to 1.
//
//   Copyright (c) 2007 Gabriel Peyre

if ~exists('options')
	options.null = 0;
end

// [vertex,face] = check_face_vertex(vertex,face);

nface = size(face,2);
n = max(max(face));

if myisfield(options, 'verb')
    verb = options.verb;
else
    verb = 0; // n>5000;
end

if ~exists('laptype')
    laptype = 'conformal';
end

select laptype
    case 'combinatorial'
        W = triangulation2adjacency(face);
    case 'distance'
        W = my_euclidean_distance(triangulation2adjacency(face),vertex);
        W(find(W>0)) = W(find(W>0)).^(-1);
        W = (W+W')/2; 
    case 'conformal'
        // conformal laplacian
        W = sparse(0); W(n,n) = 0; // allocate space for 0 matrix
		disp('Computing face ring (might take a while)');
        ring = compute_vertex_face_ring(face);
		disp('Computing weights (might take a while)');
        for i = 1:n
            if verb
//                progressbar(i,n);
            end
            for t = 1:length(ring(i).entries)
				b = ring(i).entries(t);
                // b is a face adjacent to a
                bf = face(:,b);
                // compute complementary vertices
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
                // angles
                alpha1 = myangle(vk-vi,vk-vj);
                beta1 = myangle(vj-vi,vj-vk);
                // add weight
                W(i,j) = W(i,j) + 1/tan( alpha1 ) + 1/tan( beta1 );
            end
        end
    else
        error('Unknown laptype.')
end


if myisfield(options, "normalize") 
if options.normalize==1
    W = diag(sum(W,2).^(-1)) * W;
end
end

endfunction


//////////////////////////////////////////////////////////////////////////////////
function beta1 = myangle(u,v);

du = sqrt( sum(u.^2) );
dv = sqrt( sum(v.^2) );
du = max(du,1e-15); dv = max(dv,1e-15);
beta1 = acos( sum(u.*v) / (du*dv) );

endfunction

//////////////////////////////////////////////////////////////////////////////////////////////
function W = my_euclidean_distance(A,vertex)

if size(vertex,1)<size(vertex,2)
    vertex = vertex';
end

[i,j] = find(sparse(A));
d = sum( (vertex(i,:) - vertex(j,:)).^2, 2);
W = sparse([i(:),j(:)],d);  

endfunction