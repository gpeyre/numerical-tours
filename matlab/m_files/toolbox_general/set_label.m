function set_label(xstr,ystr,zstr)

% set_label - set the label for a plot
%
%   set_label(xstr,ystr,zstr);
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin<1
    xstr = '';
end
if nargin<2
    ystr = '';
end
if nargin<3
    zstr = '';
end
set_graphic_sizes();
if not(isempty(xstr))
    h = xlabel(xstr); 
    set(h, 'FontSize', 20);
end
if not(isempty(ystr))
    h = ylabel(ystr);
    set(h, 'FontSize', 20);
end
if not(isempty(zstr))
    h = zlabel(zstr);
    set(h, 'FontSize', 20);
end