% Dual propagation.
options.end_points = [];
iterlist = .37*[.25 .5 .75 1]*n^2;
% extract first the geodesic paths
options.nb_iter_max = Inf;
[D,S] = perform_fast_marching(1./W, x0(:,1), options);
gamma = compute_geodesic(D,x0(:,2));
% iterations
clf;
for i=1:4
    options.nb_iter_max = iterlist(i);
    [D,S] = perform_fast_marching(1./W, x0, options);
    subplot(2,2,i);        
    hold on;
    imageplot(convert_distance_color(D,f));
    if i==4
        h = plot(gamma(2,:),gamma(1,:), 'k'); set(h, 'LineWidth', 2);
        % select extremal point
        u = interp2(1:n,1:n,D,gamma(2,:),gamma(1,:));
        [tmp,i] = max(u); q = gamma(:,i);
        h = plot(q(2,:),q(1,:), '.b'); set(h, 'MarkerSize', 25);
    end
    h = plot(x0(2,:),x0(1,:), '.r'); set(h, 'MarkerSize', 25);
end
