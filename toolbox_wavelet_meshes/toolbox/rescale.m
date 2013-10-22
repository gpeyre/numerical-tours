function y = rescale(x,a,b)

% rescale - rescale data in [a,b]
%
%   y = rescale(x,a,b);
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    a = 0;
end
if nargin<3
    b = 1;
end

m = min(x(:));
M = max(x(:));

y = (b-a) * (x-m)/(M-m) + a;

