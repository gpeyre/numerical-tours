function y = clamp(x,a,b)

% clamp - clamp a value
%
%   y = clamp(x,a,b);
%
% Default is [a,b]=[0,1].
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    a = 0;
end
if nargin<3
    b = 1;
end

y = max(x,a);
y = min(y,b);