function [e1,e2,l1,l2] = perform_tensor_decomp(T,order,b,c)

% perform_tensor_decomp - perform an eigendecomposition.
%
%   [e1,e2,l1,l2] = perform_tensor_decomp(T, order);
% order
%   T = perform_tensor_decomp(e1,e2,l1,l2);
%
%   e1(i,j,:) is the main eigenvector at location (i,j)
%       with associated largest eigenvalue 'l1(i,j)'.
%   e2(i,j,:) is the second eigenvector at location (i,j)
%       with associated smallest eigenvalue 'l2(i,j)'.
%
%   So you always have l1>=l2 (not in absolute value !).
%
%   If 'order'=='abs' then the the decomposition is done
%   so that abs(l1)>=abs(l2)
%
%   'T' must be a tensorial field (produced eg. by compute_hessian), 
%   so it should be symmetric.
%
%   See also: perform_tensor_decomp_3d.
%
%   Copyright (c) 2004 Gabriel Peyre

if nargin==4
    e1 = perform_tensor_recomp(T,order,b,c);
    return;
end

% retrieve the 4 entries of the tensor field
if ndims(T)==4 && size(T,3)==2 && size(T,3)==2
    K11 = T(:,:,1,1);
    K12 = T(:,:,1,2);
    K21 = T(:,:,2,1);
    K22 = T(:,:,2,2);
elseif ndims(T)==3 && size(T,3)==3
    K11 = T(:,:,1);
    K22 = T(:,:,2);
    K12 = T(:,:,3);
    K21 = K12;
elseif ndims(T)==2
    K11 = T(1,1);
    K12 = T(1,2);
    K21 = T(2,1);
    K22 = T(2,2);    
else
    error('T must be a tensor or a tensor field.');    
end

[n,p] = size(K11);

e1 = zeros(n,p,2);
e2 = zeros(n,p,2);
l1 = zeros(n,p);
l2 = zeros(n,p);

% trace/2
t = (K11+K22)/2;

a = K11 - t;
b = K12;

ab2 = sqrt(a.^2+b.^2);
l1 = ab2  + t;
l2 = -ab2 + t;

theta = atan2( ab2-a, b );

e1(:,:,1) = cos(theta);
e1(:,:,2) = sin(theta);
e2(:,:,1) = -sin(theta); 
e2(:,:,2) = cos(theta);

if nargin==2 & strcmp(order,'abs')
    % reorder the eigenvalue according to absolute value.
    A = abs(l1)>abs(l2);
    ee1 = prod_vf_sf(e1,A) + prod_vf_sf(e2,1-A);
    e2 = prod_vf_sf(e1,1-A) + prod_vf_sf(e2,A);
    e1 = ee1;
    ll1 = l1.*A + l2.*(1-A);
    l2 = l1.*(1-A) + l2.*A;
    l1 = ll1;
end



function T = perform_tensor_recomp(e1,e2,l1,l2)

% perform_tensor_recomp - create the tensor field corresponding to the given eigendecomposition.
%
%   T = perform_tensor_recomp(e1,e2,l1,l2);
%
%   'e1(i,j,:)' is the main eigenvector at location (i,j)
%       with associated largest eigenvalue 'l1(i,j)'.
%   'e2(i,j,:)' is the second eigenvector at location (i,j)
%       with associated smallest eigenvalue 'l2(i,j)'.
%
%   You have 
%       T = l1*e1*e1' + l2*e2*e2'
%
%   Copyright (c) 2004 Gabriel Peyré


T = zeros( [size(l1),2,2] );
T(:,:,1,1) = l1.*e1(:,:,1).^2 + l2.*e2(:,:,1).^2;
T(:,:,1,2) = l1.*e1(:,:,1).*e1(:,:,2) + l2.*e2(:,:,1).*e2(:,:,2);
T(:,:,2,1) = T(:,:,1,2);
T(:,:,2,2) = l1.*e1(:,:,2).^2 + l2.*e2(:,:,2).^2;
