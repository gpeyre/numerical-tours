function C = cell_add(A,B,a,b)

% cell_add - add two cell arrays
%
%   C = cell_add(A,B,a,b);
%
%   C{i} = a*A{i} + b*B{i};
%
%   Copyright (c) 2008 Gabriel Peyre


if nargin<3
    a = 1;
end
if nargin<4
    b = 1;
end

if iscell(A)
    if length(A)~=length(B)
        error('A and B must be of the same size');
    end
    for i=1:length(A)
        C{i} = cell_add(A{i},B{i},a,b);
    end
    return;
end

C = a*A + b*B;