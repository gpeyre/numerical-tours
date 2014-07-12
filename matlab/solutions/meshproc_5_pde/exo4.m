f1 = f0;
f = f0;
m = 1;
clf;
for i=1:niter
    % step in time
    [f,f1] = update(f,f1);
    if mod(i,niter/4)==1
        % display
        rho = .4;
        vmax = max(abs(f));
        a = f; a(a==max(a)) = vmax;
        a(a==min(a)) = -vmax;
        options.face_vertex_color = clamp(a,-rho*vmax,rho*vmax);
        subplot(2,2,m);
        plot_mesh(X0,F, options);
        colormap(jet(256));
        m = m+1;
    end
end
