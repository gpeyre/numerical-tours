function [f1,face1] = perform_mesh_subdivision(f, face, nsub, options)

% perform_mesh_subdivision - perfrom a mesh sub-division
%
%   [face1,f1] = perform_mesh_subdivision(f, face, nsub, options);
%
%   face is a (3,nface) matrix of original face adjacency
%   face1 is the new matrix after subdivision
%   f is a (d,nvert) matrix containing the value f(:,i) of a function
%       at vertex i on the original mesh. One should have
%           nvert=max(face(:))
%       (can be multi dimensional like point position in R^3, d=3)
%   f1 is the value of the function on the subdivided mesh.
%   
%   options.sub_type is the kind of subvision applied:
%       'linear4': 1:4 tolopoligical subivision with linear interpolation
%       'linear3': 1:3 tolopoligical subivision with linear interpolation
%       'loop': 1:4 tolopoligical subivision with loop interpolation
%       'butterfly': 1:4 tolopoligical subivision with linear interpolation
%       'sqrt3': 1:3 topological subdivision with sqrt(3) interpolation
%          (dual scheme).
%       'spherical4': 1:4 tolopoligical subivision with linear
%           interpolation and projection of f on the sphere
%       'spherical3': 1:3 tolopoligical subivision with linear
%           interpolation and projection of f on the sphere
%
%   An excellent reference for mesh subdivision is
%       Subdivision for Modeling and Animation,
%       SIGGRAPH 2000 Course notes.
%       http://mrl.nyu.edu/publications/subdiv-course2000/
%
%   The sqrt(3) subdivision is explained in
%       \sqrt{3}-subdivision, Leif Kobbelt
%       Proc. of SIGGRAPH 2000
%
%   Copyright (c) 2007 Gabriel Peyré

options.null = 0;
if nargin<2
    error('Not enough arguments');
end
if nargin==2
    nsub=1;
end

sub_type = getoptions(options, 'sub_type', '1:4');
spherical = getoptions(options, 'spherical', 0);
sanity_check = getoptions(options, 'sanity_check', 1);

switch lower(sub_type)
    case 'linear3'
        interpolation = 'linear';
        topology = 3;
    case 'linear4'        
        interpolation = 'linear';
        topology = 4;
    case 'loop'
        interpolation = 'loop';
        topology = 4;
    case 'butterfly'
        interpolation = 'butterfly';
        topology = 4;
    case 'sqrt3';
        interpolation = 'sqrt3';
        topology = 3;
    case 'spherical3'
        interpolation = 'linear';
        topology = 3;
        spherical = 1;
    case 'spherical4'
        interpolation = 'linear';
        topology = 4;
        spherical = 1;
    case '1:3'
        interpolation = 'linear';
        topology = 3;
    case '1:4'
        interpolation = 'linear';
        topology = 4;
end

if nsub==0
    f1 = f;
    face1 = face;
    return;
end

if nsub>1
    % special case for multi-subdivision
    f1 = f;
    face1 = face;
    for i = 1:nsub
         [f1,face1] = perform_mesh_subdivision(f1,face1,1, options);
    end
    return;    
end


if size(f,1)>size(f,2) && sanity_check
    f=f';
end
if size(face,1)>size(face,2) && sanity_check
    face=face';
end

m = size(face,2);
n = size(f,2);

verb = getoptions(options, 'verb', n>500);
loop_weigths = getoptions(options, 'loop_weigths', 1);

if topology==3
    f1 = ( f(:,face(1,:)) + f(:,face(2,:)) + f(:,face(3,:)))/3;
    f1 = cat(2, f, f1 );
    %%%%%% 1:3 subdivision %%%%%
    switch interpolation
        case 'linear'
            face1 = cat(2, ...
                [face(1,:); face(2,:); n+(1:m)], ...
                [face(2,:); face(3,:); n+(1:m)], ...
                [face(3,:); face(1,:); n+(1:m)] );
        case 'sqrt3'
            face1 = [];
            edge = compute_edges(face);
            ne = size(edge,2);
            e2f = compute_edge_face_ring(face);
            face1 = [];
            % create faces
            for i=1:ne
                if verb
                    progressbar(i,n+ne);
                end
                v1 = edge(1,i); v2 = edge(2,i);
                F1 = e2f(v1,v2); F2 = e2f(v2,v1);
                if min(F1,F2)<0
                    % special case
                    face1(:,end+1) = [v1 v2 n+max(F1,F2)];
                else
                    face1(:,end+1) = [v1 n+F1 n+F2];
                    face1(:,end+1) = [v2 n+F2 n+F1];
                end
            end
            % move old vertices
            vring0 = compute_vertex_ring(face);
            for k=1:n
                if verb
                    progressbar(k+ne,n+ne);
                end
                m = length(vring0{k});
               	beta = (4-2*cos(2*pi/m))/(9*m);         % warren weights
                f1(:,k) = f(:,k)*(1-m*beta) + beta*sum(f(:,vring0{k}),2);
            end
            
        otherwise 
            error('Unknown scheme for 1:3 subdivision');
    end
else
    %%%%%% 1:4 subdivision %%%%%
    i = [face(1,:) face(2,:) face(3,:) face(2,:) face(3,:) face(1,:)];
    j = [face(2,:) face(3,:) face(1,:) face(1,:) face(2,:) face(3,:)];
    I = find(i<j);
    i = i(I); j = j(I);
    [tmp,I] = unique(i + 1234567*j);
    i = i(I); j = j(I);
    ne = length(i); % number of edges
    s = n+(1:ne);

    A = sparse([i;j],[j;i],[s;s],n,n);

    % first face
    v12 = full( A( face(1,:) + (face(2,:)-1)*n ) );
    v23 = full( A( face(2,:) + (face(3,:)-1)*n ) );
    v31 = full( A( face(3,:) + (face(1,:)-1)*n ) );

    face1 = [   cat(1,face(1,:),v12,v31),...
        cat(1,face(2,:),v23,v12),...
        cat(1,face(3,:),v31,v23),...
        cat(1,v12,v23,v31)   ];
    
    
    switch interpolation
        case 'linear'            
            % add new vertices at the edges center
            f1 = [f, (f(:,i)+f(:,j))/2 ];
            
        case 'butterfly'
            
            global vring e2f fring facej;
            vring = compute_vertex_ring(face1);
            e2f = compute_edge_face_ring(face);
            fring = compute_face_ring(face);
            facej = face;
            f1 = zeros(size(f,1),n+ne);
            f1(:,1:n) = f;
            for k=n+1:n+ne
                if verb
                    progressbar(k-n,ne);
                end
                [e,v,g] = compute_butterfly_neighbors(k, n);
                f1(:,k) = 1/2*sum(f(:,e),2) + 1/8*sum(f(:,v),2) - 1/16*sum(f(:,g),2);
            end

        case 'loop'
            
            global vring e2f fring facej;
            vring = compute_vertex_ring(face1);
            vring0 = compute_vertex_ring(face);
            e2f = compute_edge_face_ring(face);
            fring = compute_face_ring(face);
            facej = face;
            f1 = zeros(size(f,1),n+ne);
            f1(:,1:n) = f;
            % move old vertices
            for k=1:n
                if verb
                    progressbar(k,n+ne);
                end
                m = length(vring0{k});
                if loop_weigths==1
                    beta = 1/m*( 5/8 - (3/8+1/4*cos(2*pi/m))^2 );   % loop original construction
                else
                    beta = 3/(8*m);         % warren weights
                end
                f1(:,k) = f(:,k)*(1-m*beta) + beta*sum(f(:,vring0{k}),2);
            end
            % move new vertices
            for k=n+1:n+ne
                if verb
                    progressbar(k,n+ne);
                end
                [e,v] = compute_butterfly_neighbors(k, n);
                f1(:,k) = 3/8*sum(f(:,e),2) + 1/8*sum(f(:,v),2);
            end
            
        otherwise 
            error('Unknown scheme for 1:3 subdivision');
    end
end

if spherical
    % project on the sphere
    d = sqrt( sum(f1.^2,1) );
    d(d<eps)=1;
    f1 = f1 ./ repmat( d, [size(f,1) 1]);
end