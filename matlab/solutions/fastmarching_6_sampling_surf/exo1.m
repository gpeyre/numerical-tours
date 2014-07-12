m = 500;
% initialize
landmarks = [100];
[D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks);
clf;
k = 1; displist = [10 50 100 m];
for i=2:m
    % select
    [tmp,landmarks(end+1)] = max(D);
    % update
    options.constraint_map = D;
    [D1,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks,options);
    D = min(D,D1);
    if i==displist(k)
        subplot(2,2,k);
        hold on;
        options.face_vertex_color = perform_hist_eq(D,'linear');
        plot_mesh(vertex,faces, options);
        colormap jet(256);
        h = plot3(vertex(1,landmarks), vertex(2,landmarks), vertex(3,landmarks), 'r.');
        set(h, 'MarkerSize', 20);
        k = k+1;
    end
end
