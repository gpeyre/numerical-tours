W = load_image('cavern');
n = size(W,1);
W = rescale(-W)+.05;

options.method = 'fm';
x0 = round([.95,.1]'*n);
D = perform_dijkstra_fm(W, x0, options);

options.order = 2;
G0 = grad(D, options);
G = G0;
% G = G0 ./ repmat( sqrt( sum(G0.^2, 3) ), [1 1 2]);

tau = 1;

x1 = [55 175]';
x1 = round([.15;.95]*n);
gamma = x1;

meth = 'nearest';
meth = 'cubic';
Geval = @(G,x)[interp2(1:n,1:n,G(:,:,1),x(2),x(1), meth); ...
             interp2(1:n,1:n,G(:,:,2),x(2),x(1), meth) ];
gamma = x1;
for i=1:15*n/tau
    gamma(:,end+1) = gamma(:,end) - tau*Geval(G, gamma(:,end));
    if norm(gamma(:,end)-x0)<1
        break;
    end
end
gamma(:,end+1) = x0;


clf; hold on;
imageplot(W); colormap gray(256);
h = plot(gamma(2,:),gamma(1,:), '.b'); set(h, 'LineWidth', 2);
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);
axis ij;