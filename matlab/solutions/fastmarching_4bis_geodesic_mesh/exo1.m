options.niter = 250*3;
options.verb = 0;
options.svg_rate = 10;
options.U = [];
[U,err,Usvg] = perform_geodesic_iterative(vertex, faces, W, I, options);
% display
t = round(linspace(1,size(Usvg,2)*.5, 4));
clf;
for kdisp=1:4
    subplot(2,2,kdisp);
    options.face_vertex_color = Usvg(:,t(kdisp));
    plot_mesh(vertex,faces,options);
    colormap jet(256);
    shading interp;
end
