function set_axis(v)

% set_axis - draw axis on/off
%
%  OFF: set_axis(0);
%  ON:  set_axis(1);
%
%   Copyright (c) 2008 Gabriel Peyre

if v==0
    axis off;
else
    axis on;
end
