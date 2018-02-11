function bound = compute_shape_boundary(M)

% compute_shape_boundary - extract boundary points of a shape
%
% bound = compute_boundary(M);
%
%   bound is the boundary of the largest connected component of the shape
%   represented by M>mean(M(:))
%
%   Copyright (c) 2009 Gabriel Peyre


M = double(M>mean(M(:)));
c = contourc(M,[.5 .5]);

b = {}; bsize = [];
while not(isempty(c))
    bsize(end+1) = c(2,1);
    b{end+1} = c(:,2:bsize(end)+1);
    c(:,1:bsize(end)+1) = [];
end
[tmp,I] = max(bsize);
bound = b{I};
bound = bound(2:-1:1,:);