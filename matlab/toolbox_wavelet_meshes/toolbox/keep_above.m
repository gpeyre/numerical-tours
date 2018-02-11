function [y,m] = keep_above(x,T)

% keep_above - keep only the coefficients above threshold T, 
%   set the rest to zero.
%
% [y,m] = keep_above(x,T);
%
%   y is the thresholded signal.
%   m is the number of coefficient above T.
%
%   Copyright (c) 2004 Gabriel Peyré

if ~iscell(x)
    % images or vectors
    I = find(abs(x)<T);
    y = x;
    y(I) = 0;
    m = sum( abs(x(:))>=T );
elseif iscell(x)
    % cell array
    for i=1:length(x)
        [y{i},m{i}] = keep_above(x{i},T);    
    end
end