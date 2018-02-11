function h = plot_sparse_diracs(x, options)

% plot_sparse_diracs - plot sparse sets of diracs
%
%   h = plot_sparse_diracs(x, options);
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin>2
    error('Works only with 2 arguments, use options.');
end

options.null = 0;
mytitle = getoptions(options, 'title', []);
delta = getoptions(options, 'delta', 1);

if iscell(x)
    clf;    
    for i=1:length(x)
        ax(i) = subplot(length(x), 1, i);
        if iscell(mytitle)
            options.title = mytitle{i};
        end
        h{i} = plot_sparse_diracs(x{i}, options);
    end
    linkaxes(ax,'x');
    return;
end

lw      = getoptions(options, 'lw', 2);
tmax  	= getoptions(options, 'tmax', .95);
val     = getoptions(options, 'val', 1e-3);
fs      = getoptions(options, 'fs', 20);


n = length(x);
xx = x; 
xx(abs(xx)>tmax)=sign(xx(abs(xx)>tmax));
xx(abs(xx)<val)=Inf;
h = stem( 0:delta:(n*delta-1), xx, 'k'); axis([0 delta*n-1 -1 1]);
set(h, 'LineWidth', lw);
set(gca, 'FontSize', fs);
if not(isempty(mytitle)) && not(iscell(mytitle))
    title(mytitle);
end
