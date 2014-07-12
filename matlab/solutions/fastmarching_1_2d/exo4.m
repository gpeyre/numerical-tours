% load
f = load_image('cavern',n);
epsilon = 1e-2;
W = epsilon + rescale(-f);
x0 = [45;280]; x1 = [275;25];
options.nb_iter_max = Inf;
options.end_points = []; % x1;
[D,S] = perform_fast_marching(1./W, x0, options);
% gradient
G0 = grad(D, options);
G = G0 ./ repmat( sqrt( sum(G0.^2, 3) ), [1 1 2]);
gamma = x1;
for i=1:1.5*n/tau
    gamma(:,end+1) = gamma(:,end) - tau*Geval(G, gamma(:,end));
    if norm(gamma(:,end)-x0)<1
        break;
    end
end
gamma(:,end+1) = x0;
% display
clf;
subplot(1,2,1);
hold on;
imageplot(repmat(f,[1 1 3]), 'Image');
h = plot(gamma(2,:),gamma(1,:), '.b'); set(h, 'LineWidth', 2);
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);
subplot(1,2,2);
hold on;
imageplot(convert_distance_color(perform_hist_eq(D,'linear'),f), 'Distance');
h = plot(gamma(2,:),gamma(1,:), '.k'); set(h, 'LineWidth', 2);
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);
