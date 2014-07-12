fj = f;
clf;
for j=1:4
    fj = reshape(fj, [length(fj)/4 4]);
    fj = mean(fj,2);
    subplot(2,2,j);
    options.face_vertex_color = fj;
    plot_mesh(vertex{end-j}, face{end-j}, options);
    view(vv);
    colormap gray(256);
    lighting none;
end
