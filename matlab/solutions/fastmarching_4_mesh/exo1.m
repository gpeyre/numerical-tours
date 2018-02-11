clf;
nblist = round( linspace(.1,1,6)*nvert );
for i = 1:length(nblist);
    options.nb_iter_max = nblist(i);
    [D,S,Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options);
    subplot(2,3,i);
    col = D; col(col==Inf) = 0;
    options.face_vertex_color = col;
    hold('on');
    plot_mesh(vertex,faces,options);
    colormap jet(256);
    % display here starting points
end
