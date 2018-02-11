function set_colormap(a)

% set_colormap - set colors for display
%
%  set_colormap('jet');
%  set_colormap('gray');
%  set_colormap('hot');
%
%  Copyright (c) 2008 Gabriel Peyre

if strcmp(a, 'jet')
    colormap jet(256); 
elseif strcmp(a, 'gray')
    colormap gray(256); 
elseif strcmp(a, 'hot')
    colormap hot(256); 
else
    error('Unknown colormap');
end