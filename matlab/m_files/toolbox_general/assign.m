function x = assign(x,I,u,  dimension)

% assign - assign values to a sub-set
%
%   assign(x,I,u)
% is equivalent to
%   x(I) = u;
% and it returns x.
%
%   assign(x,I,u, dimension)
% is equivalent to either
%   x(I,:,:) = u;   % dimension=1
%   x(:,I,:) = u;   % dimension=2
%   x(:,:,I) = u;   % dimension=3
% and it returns x.
%
%   Copyright (c) 2012 Gabriel Peyre

s = size(x);
if nargin<=3
x(I) = u;
else
    switch dimension
        case 1
            x(I,:,:) = u;
        case 2
            x(:,I,:) = u;
        case 3
            x(:,:,I) = u;
        otherwise
            error('Works only for dimension 1, 2, 3.');
    end
    
end
x = reshape(x,s);


end