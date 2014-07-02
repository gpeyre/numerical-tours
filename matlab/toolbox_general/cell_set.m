function x = cell_set(x,i,v)

% cell_set - assign value
%
%   x = cell_set(x,i,v);
% is equivalent to 
%   x{i}=v;
%
%   Copyright (c) 2008 Gabriel Peyre

x{i} = v;