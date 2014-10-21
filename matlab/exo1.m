n = 400;
landmarks = 1;
Dland = [];
k = 1;
displ = round(linspace(0,1,5)*n); displ(1) = [];
clf;
for i=1:n
    if not(isempty(Dland))
        [tmp,landmarks(end+1)] = max( min(Dland,[],2) );
    end
    [Dland(:,end+1),S,Q] = perform_fast_marching_mesh(vertex, faces, landmarks(end));
    if i==displ(k)
        options.start_points = landmarks;
        subplot(2,2,k);
        options.colorfx = 'equalize';
        plot_fast_marching_mesh(vertex,faces, min(Dland,[],2) , [], options);
        k = k+1;
    end
end
