f = f0;
k = 1; 
displist = round([.05 .1 .5 1]*niter); 
clf;
for i=1:niter
    % step
    f = f - tau*tL*f;   
    if displist(k)==i
        subplot(2,2,k);
        options.face_vertex_color = f(:);
        plot_mesh(X0,F, options);
        lighting none;
        k = k+1;
    end
end
