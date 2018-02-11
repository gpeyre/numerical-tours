function h = plot_levelset(M,t,A)

// plot_levelset - display the level set of an image
//
//   plot_levelset(M,t,A);
//
//  WARNING: NOT YET IMPLEMENTED
//
//   Copyright (c) 2008 Gabriel Peyre

if argn(2)<3
    A = M;
end

x = 1:size(M,1);
y = 1:size(M,2);

imageplot(A);
// contour(x,y,M, t, flag=[200]);
// set(h, 'LineWidth', 2);
// axis ij;


h = [];

endfunction