tlist = [.001 1 100 10000];
clf;
for k=1:length(tlist)
    t = tlist(k);
    % solve
    u = (speye(n)+t*Delta)\delta;
    % compute gradient
    g = Grad(u);
    % Normalize it to obtain
    h = -g ./ repmat( sqrt(sum(g.^2)), [3 1] );
    % Integrate it back by solving
    phi = Delta \ Div(h);
    % display
    options.face_vertex_color = DispFunc(phi);
    subplot(2,2,k);
    plot_mesh(X,F,options); 
    title(['t=' num2str(t)]);
    axis('tight');
    colormap parula(256);
end
