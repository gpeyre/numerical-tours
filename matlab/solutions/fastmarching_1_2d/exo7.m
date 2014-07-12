niter = round(linspace(.1,1,6)*n^2);
clf;
for i=1:length(niter)
    options.nb_iter_max = niter(i);
    options.end_points = [];
    [D,S] = perform_fast_marching(1./W, x0, options);    
    subplot(2,3,i);
    hold on;
    imageplot( convert_distance_color(D,f) );
    h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
end
