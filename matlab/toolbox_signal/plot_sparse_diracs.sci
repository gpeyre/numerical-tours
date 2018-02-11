function h = plot_sparse_diracs(x, options)

// plot_sparse_diracs - plot sparse sets of diracs
//
//   h = plot_sparse_diracs(x, options);
//
//   Copyright (c) 2008 Gabriel Peyre


options.null = 0;
mytitle = getoptions(options, 'title', []);
delta = getoptions(options, 'delta', 1);


lw      = getoptions(options, 'lw', 2);
tmax  	= getoptions(options, 'tmax', .95);
val     = getoptions(options, 'val', 1e-3);
fs      = getoptions(options, 'fs', 20);


n = length(x);
xx = x; 
xx(abs(xx)>tmax)=sign(xx(abs(xx)>tmax));
xx(abs(xx)<val)=0;
plot2d3( 0:delta:(n*delta-1), xx); 
axis([0 delta*n-1 -1 1]);
// set(h, 'LineWidth', lw);
// set(gca, 'FontSize', fs);
set_graphic_sizes([],fs,[]);
if not(isempty(mytitle)) & not(iscell(mytitle))
    title(mytitle);
end


endfunction