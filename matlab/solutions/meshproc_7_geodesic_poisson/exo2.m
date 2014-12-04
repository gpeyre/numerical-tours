tlist = [.001 1 100 10000];
clf;
for k=1:length(tlist)
    t = tlist(k);
    % solve
    u = (Ac+t*DeltaCot)\delta;
    % re-integrate normalized gradient
    phi = -Delta \ Div( normalize( Grad(u) ) );
    % display
    options.face_vertex_color = DispFunc(phi);
    subplot(2,2,k);
    plot_mesh(X,F,options); 
    title(['t=' num2str(t)]);
    axis('tight');
    colormap parula(256);
end
