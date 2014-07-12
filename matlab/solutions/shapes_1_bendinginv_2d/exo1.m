ms = 20;
U = geod(x);
npaths = 30;
sel = round(linspace(1,L+1,npaths+1)); sel(end) = [];
end_points = b(:,sel);
clf; hold on;
imageplot(1-S);
for i=1:npaths
    p = compute_geodesic(U,end_points(:,i));
    h = plot( p(2,:), p(1,:), 'g' ); set(h, 'LineWidth', lw);
end
h = plot(x(2),x(1), '.r'); set(h, 'MarkerSize', ms);
h = plot(end_points(2,:),end_points(1,:), '.b'); set(h, 'MarkerSize', ms);
axis ij;
