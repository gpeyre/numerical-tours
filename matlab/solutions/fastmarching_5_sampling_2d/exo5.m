% Seed random points.
p = 150;
vertex = floor(rand(2,p)*(n-1))+1;
disp_list = [1,2,5,30]; k = 1;
clf; hold on;
for i=1:max(disp_list)
    % Compute Voronoi partition.
    [D,Z,Q] = perform_fast_marching(1./W, vertex);
    if i==disp_list(k)
        subplot(2,2,k);
        imageplot(Q'); hold on;
        h = plot(vertex(1,:), vertex(2,:), 'k.');
        set(h, 'MarkerSize', 15);
        colormap(jet(256));
        k = k+1;
    end
    % Re-center each point at the barycenter of its cell.    
    for i=1:p
        [x,y] = ind2sub(size(W), find(Q==i));
        vertex(:,i) = [mean(x);mean(y)];
    end
end
