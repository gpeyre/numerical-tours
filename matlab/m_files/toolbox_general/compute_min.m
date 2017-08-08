function [Y,I] = compute_min(X,d)

% compute_min - compute min along dimension d
%
%   [Y,I] = compute_min(X,d);
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin<2
    d = 1;
end

[Y,I] = min(X,[],d);