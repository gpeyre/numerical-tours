function [Y,I] = compute_max(X,d)

% compute_max - compute maximum along dimension d
%
%   [Y,I] = compute_max(X,d);
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin<2
    d = 1;
end

[Y,I] = max(X,[],d);