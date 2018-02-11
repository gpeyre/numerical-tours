function set_label(xstr,ystr,zstr)

// set_label - set the label for a plot
//
//   set_label(xstr,ystr);
//
//   Copyright (c) 2008 Gabriel Peyre

if argn(2)<1
    xstr = '';
end
if argn(2)<2
    ystr = '';
end
if argn(2)<3
    zstr = '';
end
set_graphic_sizes();
if not(isempty(xstr))
    h = xlabel(xstr); 
end
if not(isempty(ystr))
    h = ylabel(ystr);
end
if not(isempty(zstr))
    h = zlabel(zstr);
end

endfunction