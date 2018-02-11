function f = subselectdim(f,sel,d)

% subselectdim - select along a dimension
%
%   f = subselectdim(f,sel,d);
%
% |subselectdim(A,sel,1)| is equivalent to |A(sel,:,:)|
% while |subselectdim(A,sel,2)| is equivalent to |A(:,sel,:)|.
%
% Copyright (c) 2009 Gabriel Peyre

switch d
    case 1
        f = f(sel,:,:,:,:,:,:,:);
    case 2
        f = f(:,sel,:,:,:,:,:,:);
    case 3
        f = f(:,:,sel,:,:,:,:,:);
    case 4
        f = f(:,:,:,sel,:,:,:,:);
    case 5
        f = f(:,:,:,:,sel,:,:,:);
    case 6
        f = f(:,:,:,:,:,sel,:,:);
    case 7
        f = f(:,:,:,:,:,:,sel,:);
    case 8
        f = f(:,:,:,:,:,:,:,sel);
    otherwise
        error('Not implemented');
end
