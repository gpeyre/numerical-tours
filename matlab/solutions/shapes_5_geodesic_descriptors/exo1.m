p = sum(M(:)==1);
qlist = ceil([.1 .3 .7 1]*p);
clf;
for i=1:4
    options.nb_iter_max = qlist(i);
    D = perform_fast_marching(W, start_points, options);
    subplot(2,2,i);
    display_shape_function(perform_hist_eq(D, 'linear')); 
    hold on; 
    h = plot(bound(2,:), bound(1,:), 'k');
    set(h, 'LineWidth', 2);
    axis('ij');
end
options.nb_iter_max = Inf;
