clf;
options.W = W;
nblist = round( [.05 .15 .4 .6]*nvert );
for i = 1:length(nblist);
    options.nb_iter_max = nblist(i);
    [D,S,Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options);
    subplot(2,2,i);
    plot_fast_marching_mesh(vertex,faces, D, [], options);
end
