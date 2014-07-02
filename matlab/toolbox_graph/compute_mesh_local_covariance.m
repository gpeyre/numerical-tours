function [C,U,D] = compute_mesh_local_covariance(vertex,face,f,options)

% compute_mesh_local_covariance - compute the local covariance of a vector field
%
%   [C,U,D] = compute_mesh_local_covariance(vertex,face,f,options);
%
%   f(:,i) is the vector value at some vertex i.
%   C(:,:,i) is the covariance of f around vertex i.
%
%   The number of vertices over which the covariance is computed
%   is determined by options.covariance_smoothing.
%
%   U(:,1) is an approximate of the normal.
%   U(:,2) is an approximate of the 1st tangent to the surface.
%   U(:,3) is an approximate of the 1st tangent to the surface.
%
%   References for the methods include:
%   Integral Invariants for Robust Geometry Processing 
%   Helmut Pottmann, Johannes Wallner, Qixing Huang, and Yong-Liang Yang 
% and
%   Robust Feature Detection and Local Classification for Surfaces Based on Moment Analysis
%   Ulrich Clarenz, Martin Rumpf, Alexandru Telea 
%
%   Copyright (c) 2007 Gabriel Peyre


[vertex,face] = check_face_vertex(vertex,face);
if size(f,1)>size(f,2)
    f = f';
end

n = size(vertex,2);
m = size(face,2);
d = size(f,1);

options.niter_averaging = getoptions(options, 'covariance_smoothing', 10);

C = zeros(d,d,n);
for i=1:d
    for j=1:i
        ai = perform_mesh_smoothing(face,vertex,f(i,:)',options);
        aj = perform_mesh_smoothing(face,vertex,f(j,:)',options);
        aij = perform_mesh_smoothing(face,vertex, f(i,:)' .* f(j,:)',options);
        C(i,j,:) =  reshape(aij,[1 1 n]) - reshape(ai,[1 1 n]).*reshape(aj,[1 1 n]);
        C(j,i,:) = C(i,j,:);
    end
end

if nargout>1
    % extract eigenvectors and eigenvalues
    U = zeros(d,d,n);
    D = zeros(d,n);
    for k=1:n
        progressbar(k,n);
        [u,d] = eig(C(:,:,k));
        d = real(diag(d));
        % sort acording to [norma,min curv, max curv]
        [tmp,I] = sort(abs(d));
        D(:,k) = d(I);
        U(:,:,k) = real(u(:,I));
    end

    % try to re-orient the normals
    normal = compute_normal(vertex,face);
    Normal = squeeze(U(:,1,:));
    s = sign( sum(Normal.*normal,1) );
    Normal = Normal .* repmat(s, 3,1);
    U(:,1,:) = reshape(Normal, [3 1 n]);
end
    