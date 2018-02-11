clf;
klist = [1 2 4 8];
i = 1;
f1 = f;
for k=1:max(klist)
    f1 = tW*f1;
    if k==klist(i)
        options.face_vertex_color = f1(:);
        subplot(2,2,i);
        plot_mesh(X0,F, options);
        lighting none;
        i = i+1;
    end
end
options.face_vertex_color = [];
