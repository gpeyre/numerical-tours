function set_colormap(a)

// set_colormap - set colors for display
//
//  set_colormap('jet');
//  set_colormap('gray');
//  set_colormap('hot');
//
//  Copyright (c) 2008 Gabriel Peyre

f=gcf(); 
if strcmp(a, 'jet')
    f.color_map=jetcolormap(256); 
elseif strcmp(a, 'gray')
    f.color_map=graycolormap(256); 
elseif strcmp(a, 'hot')
    f.color_map=hotcolormap(256); 
else
    error('Unknown colormap');
end

endfunction