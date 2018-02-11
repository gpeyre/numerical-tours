function set_graphic_sizes(h,fs,lw)

% set_graphic_sizes - enlarge the size of the fonts
%
%   set_graphic_sizes(h,fs,lw);
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin==0
    h = [];
end
if nargin<2
    fs = 20;
end
if nargin<3
    lw = 2;
end

if not(isempty(h))
    set(h, 'LineWidth', lw);
end
set(gca, 'FontSize', fs);