clf;
for i=1:4
    subplot(2,2,i);
    plot_triangulation(vertex_svg{i},faces_svg{i}, M);
end
colormap jet(256);
