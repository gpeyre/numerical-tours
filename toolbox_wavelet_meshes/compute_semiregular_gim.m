function [vertex,face,vertex0] = compute_semiregular_gim(M,J,options)

% compute_semiregular_gim - compute a semi-regular mesh from a GIM
%
%   [vertex,face,vertex0] = compute_semiregular_gim(M,J,options);
%
%   M is an (n,n,3) geometry image.
%   vertex and face is a cell array of semi regular meshes.
%
%   Copyright (c) 2007 Gabriel Peyre

vertex0 = {}; face = {};
[vertex0{1}, face{1}] = compute_base_mesh('oct');
options.base_mesh = 'oct';
for j=2:J
    [vertex0{j},face{j}] = perform_mesh_subdivision(vertex0{j-1},face{j-1},1);
end


for j=1:J
    pos = vertex0{j};
    uv = pos(1:2,:);
    uv1 = uv(:,pos(3,:)>0);
    % un-wrapping #1
    I = find(uv1(1,:)>=0 & uv1(2,:)>=0);
    uv1(:,I) = 1 - uv1(2:-1:1,I);
    % un-wrapping #2
    I = find(uv1(1,:)>=0 & uv1(2,:)<=0);
    uv1(:,I) = uv1(2:-1:1,I);
    uv1(1,I) = uv1(1,I)+1;
    uv1(2,I) = uv1(2,I)-1;
    % un-wrapping #3
    I = find(uv1(1,:)<=0 & uv1(2,:)<=0);
    uv1(:,I) = - 1 - uv1(2:-1:1,I);
    % un-wrapping #4
    I = find(uv1(1,:)<=0 & uv1(2,:)>=0);
    uv1(:,I) = uv1(2:-1:1,I);
    uv1(1,I) = uv1(1,I)-1;
    uv1(2,I) = uv1(2,I)+1;
    % assign
    uv(:,pos(3,:)>0) = uv1;
    % rescale
    n = size(uv,2); % number of grid points
    q = size(M,1);
    uv = (uv+1)/2; % uv = uv*(1-1/q);
    % perform interpolation
    x = linspace(0,1,q);
    vertex{j} = zeros(3,n);
    for i=1:3
        vertex{j}(i,:) = interp2(x,x,M(:,:,i), uv(1,:), uv(2,:));
    end
end