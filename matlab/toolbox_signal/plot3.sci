function h = plot3(x,y,z,s)

// plot3 - 3D plot
//
//  plot3(x,y,z,s);
//
//  With s='k' or 'k.' for points.
//
//  Copyright (c) 2008 Gabriel Peyre

if argn(2)<4
    s = 'k';
end

x = x(:);
y = y(:);
z = z(:);
if strcmp(s, 'k.') | strcmp(s, 'b.') | strcmp(s, '.')
    sty = -4;
else
    sty = 1;
end

param3d1(x,y, list(z(:), sty));

h = gce();

endfunction