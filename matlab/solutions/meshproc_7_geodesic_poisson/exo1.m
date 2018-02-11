t = .1;
nlist = [1 3 10 50]; k = 1;
clf;
u = delta;
for i=1:max(nlist)
    % solve
    u = (Ac+t*DeltaCot)\u;
    % display
    if k<=4 && i==nlist(k)
        subplot(2,2,k); hold on;
        options.face_vertex_color = u;
        plot_mesh(X,F,options);
        axis('tight');
        colormap parula(256);
        title(['Iter ' num2str(i)]);
        drawnow;
        k = k+1;
    end
end
