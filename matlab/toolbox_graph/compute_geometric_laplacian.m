function L = compute_geometric_laplacian(vertex,face,type)

% compute_geometric_laplacian - return a laplacian
%   of a given triangulation (can be combinatorial or geometric).
%
%   L = compute_geometric_laplacian(vertex,face,type);
%
%   Type is either : 
%       - 'combinatorial' : combinatorial laplacian, doesn't take into
%       acount geometry.
%       - 'conformal' : conformal (i.e. harmonic) weights, 
%       - 'authalic' : area-preserving weights (not implemented yet).  
%       - 'distance' : L(i,j) = 1/|xi-xj|^2
%
%       Reference: 
%           M.S.Floater and K.Hormann
%           Recent Advances in Surface Parameterization
%           Multiresolution in geometric modelling
%           <http://vcg.isti.cnr.it/~hormann/papers/survey.pdf>
%
%   Copyright (c) 2003 Gabriel Peyre


error('Not used anymore');

[vertex,face] = check_face_vertex(vertex,face);

nface = size(face,1);
n = max(max(face));

if nargin<3
    type = 'conformal';
end

if strcmp(lower(type),'combinatorial')
    L = compute_combinatorial_laplacian( triangulation2adjacency(face) );
    return;
end
if strcmp(lower(type),'distance')
    D = build_euclidean_weight_matrix(triangulation2adjacency(face),vertex,0);
    D(D>0) = 1/D(D>0).^2;
    L = diag( sum(D) ) - D;
    L = (L+L')/2; 
    return;
end
if strcmp(lower(type),'authalic')
    error('Not implemented');
end
if not(strcmp(lower(type),'conformal'))
    error('Unknown laplacian type.');
end

% conformal laplacian
L = sparse(n,n);

ring = compute_vertex_face_ring(face);

for i = 1:n
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
        L(i,j) = L(i,j) + cot( alpha );
        L(i,k) = L(i,k) + cot( beta );
    end
end

L = L - diag(sum(L,2));


return;

%% old code

if strcmp(lower(type),'combinatorial')
    L = compute_combinatorial_laplacian( triangulation2adjacency(face) );
    
elseif strcmp(lower(type),'conformal') || strcmp(lower(type),'authalic')
    if nargin<4
        disp('--> Computing 1-ring.');
        ring = compute_vertex_ring( face );
    end
    disp('--> Computing laplacian.');
    for i=1:n
        vi = vertex(i,:);
        r = ring{i};
        if r(end)==-1
            % no circularity
            s = length(r)-1;
            r = [r(1), r(1:(end-1)), r(end-1)];
        else
            % circularity
            s = length(r);
            r = [r(end), r, r(1)];
        end
        % circulate on the 1-ring        
        for x = 2:(s+1)        
            j = r(x);            
            if L(i,j)==0            
                gche = r(x-1);                
                drte = r(x+1);                
                vj = vertex(j,:);                
                v1 = vertex(gche,:);                
                v2 = vertex(drte,:);                
                % we use cot(acos(x))=x/sqrt(1-x^2)
                if strcmp(lower(type),'conformal')
                    d1 = sqrt(dot(vi-v2,vi-v2));                
                    d2 = sqrt(dot(vj-v2,vj-v2));                
                    if d1>eps && d2>eps                
                        z = dot(vi-v2,vj-v2)/( d1*d2 );                    
                        L(i,j) = L(i,j) + z/sqrt(1-z^2);                    
                    end                
                    d1 = sqrt(dot(vi-v1,vi-v1));                
                    d2 = sqrt(dot(vj-v1,vj-v1));                
                    if d1>eps && d2>eps                
                        z = dot(vi-v1,vj-v1)/( d1*d2 );                    
                        L(i,j) = L(i,j) + z/sqrt(1-z^2);                    
                    end         
                else
                    d1 = sqrt(dot(vi-vj,vi-vj));                
                    d2 = sqrt(dot(v2-vj,v2-vj));                
                    if d1>eps && d2>eps                
                        z = dot(vi-vj,v2-vj)/( d1*d2 );                    
                        L(i,j) = L(i,j) + z/sqrt(1-z^2);                    
                    end
                    d1 = sqrt(dot(vi-vj,vi-vj));
                    d2 = sqrt(dot(v1-vj,v1-vj));                 
                    if d1>eps && d2>eps                
                        z = dot(vi-vj,v1-vj)/( d1*d2 );                    
                        L(i,j) = L(i,j) + z/sqrt(1-z^2);                    
                    end
                    if d1>eps
                        L(i,j) = L(i,j) / (d1*d1);
                    end
                end
                if 0    % uncomment for symmeterization
                if L(j,i)==0
                    L(j,i) =  L(i,j);
                else
                    L(j,i) =  (L(j,i)+L(i,j))/2;
                end
                end
            end            
        end        
    end
            
    for i=1:n    
        L(i,i) = -sum( L(i,:) );        
    end
else
    error('Unknown type.');        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function beta = myangle(u,v);

du = sqrt( sum(u.^2) );
dv = sqrt( sum(v.^2) );
du = max(du,eps); dv = max(dv,eps);
beta = acos( sum(u.*v) / (du*dv) );