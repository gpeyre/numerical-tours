function y = perform_radon_transform(x, tlist, direction, rotation)

% perform_radon_transform - Radon transform
%
%   y = perform_radon_transform(x, tlist, direction, rotation);
%
%   tlist is a vector of m angles.
%   If direction==1, compute the Radon transform.
%       x is of size (n,n) and y of size (n,m)
%   If direction==-1, compute the adjoint of Radon transform (not the inverse !).
%       x is of size (n,m) and y of size (n,n)
%
%   Copyright (c) 2011 Gabriel Peyre

n = size(x,1);
m = length(tlist);

if nargin<3
    t = linspace(-1,1,n);
    [Y,X] = meshgrid(t,t);
    rotation = @(f,t)interp2(Y,X,f, sin(t)*X + cos(t)*Y, cos(t)*X - sin(t)*Y, 'cubic', 0); 
end

if direction==1
    %% Radon %%
    y = zeros(n,m);
    for i=1:m
        y(:,i) = sum(rotation(x, tlist(i)), 2);
    end
else
    %% Adjoint %%
    y = zeros(n);
    for i=1:m
        y = y + rotation(repmat(x(:,i), [1 n]), -tlist(i));
    end
end