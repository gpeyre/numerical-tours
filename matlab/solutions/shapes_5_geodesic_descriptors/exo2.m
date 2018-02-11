npaths = 30;
sel = round(linspace(1,nbound+1,npaths+1)); sel(end) = [];
end_points = bound(:,sel);
clf; hold on;
imageplot(1-M);
for i=1:npaths
    p = compute_geodesic(D,end_points(:,i));
    h = plot( p(2,:), p(1,:), 'g' ); set(h, 'LineWidth', lw);
end
h = plot(start_points(2),start_points(1), '.r'); set(h, 'MarkerSize', ms);
h = plot(end_points(2,:),end_points(1,:), '.b'); set(h, 'MarkerSize', ms);
axis ij;
