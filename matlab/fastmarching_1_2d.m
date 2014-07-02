%% Fast Marching in 2D
% This tour explores the use of Fast Marching methods in 2-D.

perform_toolbox_installation('signal', 'general', 'graph');

%CMT
rep = 'results/fastmarching_2d/';
if not(exist(rep))
    mkdir(rep);
end
%CMT

%% Shortest Path for Isotropic Metrics
% Shortest paths are 2D curves that minimize a weighted length according to
% a given metric \(W(x)\) for \(x \in [0,1]^2\).
% The metric is usually computed from an input image \(f(x)\).

%%
% The length of a curve \( t \in [0,1] \mapsto \gamma(t) \in [0,1]^2 \) is 
% \[ L(\gamma) = \int_0^1 W(\gamma(t)) \norm{\gamma'(t)} \text{d} t. \]

%%
% Note that \(L(\gamma)\) is invariant under re-parameterization of the
% curve \(\gamma\).

%%
% A geodesic curve \(\gamma\) between two points \(x_0\) and \(x_1\) has minimum
% length among curves joining \(x_0\) and \(x_1\),
% \[ \umin{\ga(0)=x_0, \ga(1)=x_1} L(\ga). \]
% A shortest curve thus tends to pass in areas where \(W\) is small.


%% 
% The geodesic distance between the two points is then 
% \(d(x_0,x_1)=L(\gamma)\) is the geodesic distance according to the metric \(W\).

%% Pixel values-based Geodesic Metric
% The geodesic distance map \(D(x)=d(x_0,x)\) to a fixed starting point \(x_0\)
% is the unique viscosity solution of 
% the Eikonal equation 
% \[ \norm{ \nabla D(x)} = W(x) \qandq D(x_0)=0. \]

%%
% This equation can be solved numerically in \(O(N \log(N))\) operation on a discrete
% grid of \(N\) points.


%%
% We load the input image \(f\).

clear options;
n = 300;
name = 'road2';
f = rescale( load_image(name, n) );

%%
% Display the image.

clf;
imageplot(f);

%%
% Define start and end points \(x_0\) and \(x_1\) (note that you can use your own points).

x0 = [14;161];
x1 = [293;148];

%%
% The metric is defined according to \(f\) in order to be low at pixel
% whose value is close to \(f(x)\). A typical example is 
% \[ W(x) = \epsilon + \abs{f(x_0)-f(x)} \]
% where the value of \( \epsilon>0 \) should be increased in order to
% obtain smoother paths.

epsilon = 1e-2;
W = epsilon + abs(f-f(x0(1),x0(2)));

%%
% Display the metric \(W\).

clf;
imageplot(W);

%CMT
imwrite(rescale(f), [rep 'road-image.png'], 'png');
imwrite(rescale(W), [rep 'road-metric.png'], 'png');
%CMT

%%
% Set options for the propagation: infinite number of iterations, and stop
% when the front hits the end point.

options.nb_iter_max = Inf;
options.end_points = x1;

%%
% Perform the propagation, so that \(D(a,b)\) is the geodesic distance
% between the pixel \(x_1=(a,b)\) and the starting point \(x_0\).
% Note that the function |perform_fast_marching| takes as input the inverse
% of the metric \(1/W(x)\).

[D,S] = perform_fast_marching(1./W, x0, options);

%% 
% Display the propagated distance map \(D\).
% We display in color the distance map in areas where the front has
% propagated, and leave in black and white the area where the front did not
% propagate.

clf;
hold on;
imageplot( convert_distance_color(D,f) );
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);


%EXO
%% Using |options.nb_iter_max|, display the progressive propagation.
%% This corresponds to displaying the front 
%% \( \enscond{x}{D(x) \leq t} \) for various arrival times \(t\).
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
%EXO

%% Geodesic Curve Extraction
% Once the geodesic distance map \(D(x)\) to a starting point \(x_0\) is
% computed, the geodesic curve between any point \(x_1\) and \(x_0\)
% extracted through gradient descent
% \[ \ga'(t) = - \eta_t \nabla D(\ga(t)), \]
% where \(\eta_t>0\) controls the parameterization speed of the resulting
% curve. To obtain unit speed parameterization, one can use \(\eta_t =
% \norm{\nabla D(\ga(t))}^{-1}\).

%%
% Recompute the geodesic distance map \(D\) on the whole grid.

options.nb_iter_max = Inf;
options.end_points = [];
[D,S] = perform_fast_marching(1./W, x0, options);

%%
% Display \(D\).

clf;
imageplot(D);
colormap jet(256);

%%
% Compute the gradient \(G_0(x) = \nabla D(x) \in \RR^2\) of the distance map. Use centered differences.

options.order = 2;
G0 = grad(D, options);

%%
% Normalize the gradient to obtained \(G(x) = G_0(x)/\norm{G_0(x)}\), in order to have unit speed geodesic curve (parameterized
% by arc length).

G = G0 ./ repmat( sqrt( sum(G0.^2, 3) ), [1 1 2]);

%%
% Display \(G\).

clf;
imageplot(G);
colormap jet(256);

%%
% The geodesic is then numerically computed using a discretized gradient
% descent, which defines a discret curve \( (\ga_k)_k \) using 
% \[ \ga_{k+1} = \ga_k - \tau G(\ga_k) \]
% where \(\ga_k \in \RR^2\) is an approximation of \(\ga(t)\) at time
% \(t=k\tau\), and the step size \(\tau>0\) should be small enough.

%%
% Step size \(\tau\) for the gradient descent.

tau = .8;

%%
% Initialize the path with the ending point.

gamma = x1;

%%
% Define a shortcut to interpolate \(G\) at a 2-D points.
% _Warning:_ the |interp2| switches the role of the axis ...

Geval = @(G,x)[interp2(1:n,1:n,G(:,:,1),x(2),x(1)); ...
             interp2(1:n,1:n,G(:,:,2),x(2),x(1)) ];

%%
% Compute the gradient at the last point in the path, using interpolation.

g = Geval(G, gamma(:,end));

%% 
% Perform the descent and add the new point to the path.

gamma(:,end+1) = gamma(:,end) - tau*g;

%EXO
%% Perform the full geodesic path extraction by iterating the gradient
%% descent. You must be very careful when the path become close to
%% \(x_0\), because the distance function is not differentiable at this
%% point. You must stop the iteration when the path is close to \(x_0\).
gamma = x1;
for i=1:1.5*n/tau
    gamma(:,end+1) = gamma(:,end) - tau*Geval(G, gamma(:,end));
    if norm(gamma(:,end)-x0)<1
        break;
    end
end
gamma(:,end+1) = x0;
%EXO

%%
% Display the curve on the image background.

clf; hold on;
imageplot(f);
h = plot(gamma(2,:),gamma(1,:), '.b'); set(h, 'LineWidth', 2);
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);
axis ij; 


%%
% Display the curve on the distance background.

clf; hold on;
imageplot(D); colormap jet(256);
h = plot(gamma(2,:),gamma(1,:), '.b'); set(h, 'LineWidth', 2);
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);
axis ij;

%EXO
%% Study the influence of the \(\epsilon\) parameter.
elist = [.1 .5 1 10];
clf;
for ie=1:length(elist)
    epsilon = elist(ie);
    W = epsilon + abs(f-f(x0(1),x0(2)));
    [D,S] = perform_fast_marching(1./W, x0, options);
    G0 = grad(D, options);
    G = G0 ./ repmat( sqrt( sum(G0.^2, 3) ), [1 1 2]);
    % 
    gamma = x1;
    for i=1:1.5*n/tau
        gamma(:,end+1) = gamma(:,end) - tau*Geval(G, gamma(:,end));
        if norm(gamma(:,end)-x0)<1
            break;
        end
    end
    gamma(:,end+1) = x0;
    %
    subplot(2,2,ie); hold on;
    imageplot(f);
    h = plot(gamma(2,:),gamma(1,:), '.b'); set(h, 'LineWidth', 2);
    h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
    h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);
    axis ij; 
    title(['\epsilon=' num2str(epsilon)]);
end
%EXO

%EXO
%% Perform the shortest path 
%% extraction for various images such as 'cavern' or 'mountain'.
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
%EXO


%% Edge-based Geodesic Metric
% It is possible to extract the boundary of an object using shortest paths
% that follows region of high gradient.

%%
% First we load an image \(f\).

n = 256;
name = 'cortex';
f = rescale( sum(load_image(name,n),3) );

%%
% Display it.

clf;
imageplot(f);

%%
% An edge-attracting potential \(W(x)\) should be small 
% in regions of high gradient. A popular choice is
% \[ W(x) = \frac{1}{\epsilon + G_\si \star G(x)}
%  \qwhereq G(x) = \norm{\nabla f(x)}, \]
% and where \(G_\si\) is a Gaussian kernel of variance \(\si^2\).

%%
% Compute the gradient norm \(G(x)\).

G = grad(f,options);
G = sqrt( sum(G.^2,3) );

%%
% Smooth it by \(G_\si\).

sigma = 3;
Gh = perform_blurring(G,sigma);

%%
% Display the smoothed gradient \( G \star G_\si \).

clf;
imageplot(Gh);


%%
% Compute the metric.

epsilon = 0.01;
W = 1./( epsilon + Gh );


%%
% Display it.

clf;
imageplot(W);

%%
% Set two starting point \( \Ss = \{x_0^1,x_0^2\} \) (you can use other points).

x0 = [ [136;53] [123;205]];

%%
% Compute the Fast Marching from these two base points.

options.nb_iter_max = Inf;
options.end_points = [];
[D,S,Q] = perform_fast_marching(1./W, x0, options);

%%
% Display the geodesic distance (with color normalization).

clf; hold on;
imageplot( perform_hist_eq(D,'linear') );
h = plot(x0(2,:),x0(1,:), '.r'); set(h, 'MarkerSize', 25);
colormap jet(256);

%%
% The Voronoi segmentation associated to \(\Ss\) is 
% \[ \Cc_i = \enscond{x}{ \forall j \neq i, \; d(x_0^i,x) \leq d(x_0^j,x) }. \]

%%
% This Voronoi segmentation is computed during the Fast Marching
% propagation and is encoded in the partition function \(Q(x)\)
% using \(\Cc_i = \enscond{x}{Q(x)=i}\).

%%
% Display the distance and the Voronoi segmentation.

clf; hold on;
A = zeros(n,n,3); A(:,:,1) = rescale(Q); A(:,:,3) = f;
imageplot(A);
h = plot(x0(2,:),x0(1,:), '.g'); set(h, 'MarkerSize', 25);

%EXO
%% Extract the set of points that are along the boundary of the Voronoi
%% region. This corresponds for instance to the points of the region 
%% \( \enscond{x}{Q(x)=1} \)
%% that have one neighbor inside the region
%% \( \enscond{x}{Q(x)=2} \). 
%% Compute the geodesic distance \(D(x)\) at these points, and choose two points
%% \(a\) and \(b\) on this boundary that have small values of \(D\).
% Hint: you can use a convolution |U=conv2(double(Q==2),h,'same')| with a
% well chose kernel |h| to located the points |U>0| with at least 1
% neighbor.
h = [0 1 0; 1 0 1; 0 1 0];
B = (Q==1) & (conv2(double(Q==2),h,'same')>0);
U = find(B);
[xa,xb] = ind2sub(size(f),U);
[xa,I] = sort(xa); xb = xb(I); U = U(I);
dU = D(U);
k = [65 259];
x1 =[ [xa(k(1));xb(k(1))] [xa(k(2));xb(k(2))] ];
%
clf;
% subplot(2,1,1);
hold on;
imageplot(A, 'Boundary points');
h = plot(x0(2,:),x0(1,:), '.g'); set(h, 'MarkerSize', 25);
h = plot(xb,xa, 'g'); set(h, 'LineWidth', 2);
h = plot(x1(2,:), x1(1,:), '.b'); set(h, 'MarkerSize', 25);
%
if 0
subplot(2,1,2);
hold on;
h = plot( [k(1) k(1)], [min(dU) max(dU)], 'r:' ); set(h, 'LineWidth', 2);
h = plot( [k(2) k(2)], [min(dU) max(dU)], 'r:' ); set(h, 'LineWidth', 2);
h = plot(dU); axis('tight');title('D along the boundary'); set(h, 'LineWidth', 2);
end
%EXO

%EXO
%% Extract the geodesics joining \(a\) and \(b\) to the two starting points
%% (this makes 4 geodesic curves). Use them to perform segmentation.
options.nb_iter_max = Inf;
options.end_points = [];
tau = .8;
curves = {};
for i=1:2
    % FM
    [D,S,Q] = perform_fast_marching(1./W, x0(:,i), options);
    % gradient
%    D1 = D; D1(D1==Inf) = max(D1(D1~=Inf));
    G = grad(D, options);
    G = G ./ repmat( sqrt( sum(G.^2, 3)+1e-9 ), [1 1 2]);
    % Geodesic
    for j=1:2
        gamma = x1(:,j);
        for iter=1:1.5*n/tau
            g = [interp2(1:n,1:n,G(:,:,1),gamma(2,end),gamma(1,end)); ...
                interp2(1:n,1:n,G(:,:,2),gamma(2,end),gamma(1,end)) ];
            gamma(:,end+1) = gamma(:,end) - tau*g;
            if norm(gamma(:,end)-x0(:,i))<1
                break;
            end
        end
        gamma(:,end+1) = x0(:,i);
        curves{end+1} = gamma;
    end    
end
% display the curves
col = {'r', 'g', 'b', 'c'};
clf;
hold on;
imageplot(A, 'Boundary points');
for i=1:length(curves)
    c = curves{i};
    h = plot(c(2,:),c(1,:), col{i}); set(h, 'LineWidth', 2);    
end
h = plot(x0(2,:),x0(1,:), '.g'); set(h, 'MarkerSize', 25);
h = plot(x1(2,:), x1(1,:), '.b'); set(h, 'MarkerSize', 25);
%EXO


%% Vessel Segmentation and Centerline Extraction
% One can extract a network of geodesic curve starting from a central point
% to detect vessels in medical images.

%%
% Load an image. This image is extracted from the 
%  <http://www.isi.uu.nl/Research/Databases/DRIVE/ DRIVE database> of
%  retinal vessels.

n = 256;
name = 'vessels';
f = rescale(load_image(name, n));

%%
% Display it.

clf;
imageplot(f);

%%
% We clean the image by substracting the smoothly varying background
% \[ f_1 = f - G_\si \star f, \]
% where \(G_\si\) is a Gaussian kernel of variance \(\si^2\).
% Computing \(f_1\) corresponds to a high pass filtering.

sigma = 20;
f1 = perform_blurring(f,sigma) - f;

%%
% Display this normalized image.

clf;
imageplot(f1);

%%
% We compute a metric tthat is small for large values of \(f_1\):
% \[ W(x) = \epsilon +  \abs{f_1(x)-c} 
%       \qwhereq c = \umax{x} f_1(x). \]

c = max(f1(:));
epsilon = 1e-2;
W = epsilon + abs(f1-c);

%%
% Display the metric.

clf,
imageplot(W);

%%
% Select a central point \(x_0\) for the network.

x0 = [142;226];

%EXO
%% Perform partial propagations from \(x_0\).
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
%EXO

%EXO
%% Extract geodesics joining several points \(x_1\) to the central point
%% \(x_0\).
x1 = [ [175;5] [21;42] [48;133] [244;78] [191;40] ...
         [100;13] [66;42] [183;66] [220;117]];
% gradient
D1 = D; D1(D1==Inf) = max(D1(D1~=Inf));
G = grad(D1, options);
G = G ./ repmat( sqrt( sum(G.^2, 3)+1e-9 ), [1 1 2]);
% extract centerlines
curves = {};
for k=1:size(x1,2)
    % extract curve
    gamma = x1(:,k);
    for iter=1:1.5*n/tau
        g = [interp2(1:n,1:n,G(:,:,1),gamma(2,end),gamma(1,end)); ...
            interp2(1:n,1:n,G(:,:,2),gamma(2,end),gamma(1,end)) ];
        gamma(:,end+1) = clamp( gamma(:,end) - tau*g, 1,n );
        if norm(gamma(:,end)-x0)<1
            break;
        end
    end
    gamma(:,end+1) = x0;
    curves{end+1} = gamma;
end
% display the curves
clf;
hold on;
imageplot(f, 'Boundary points');
for i=1:length(curves)
    c = curves{i};
    h = plot(c(2,:),c(1,:), 'r'); set(h, 'LineWidth', 2);    
end
h = plot(x0(2,:),x0(1,:), '.g'); set(h, 'MarkerSize', 25);
h = plot(x1(2,:), x1(1,:), '.b'); set(h, 'MarkerSize', 25);
%EXO



%% Dual Propagation
% In order to speed up geodesic extraction, one can perform the propagation
% from both the start point \(x_0^1\) and end point \(x_0^2\).

%%
% Boundary points.

x0 = [[143;249] [174;9]];

%EXO
%% Perform the dual propagation, and stop it when the front meet.
%% Extract the two half geodesic curves.
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
%EXO


%CMT
iterlist = .37*linspace(0,1,6)*n^2; iterlist(1) = [];
for i=1:length(iterlist)
    options.nb_iter_max = iterlist(i);
    [D,S] = perform_fast_marching(1./W, x0, options);
    clf; hold on;
    imageplot(convert_distance_color(D,f));
    if i==length(iterlist)
        h = plot(gamma(2,:),gamma(1,:), 'k'); set(h, 'LineWidth', 2);
        % select extremal point
        u = interp2(1:n,1:n,D,gamma(2,:),gamma(1,:));
        [tmp,v] = max(u); q = gamma(:,v);
        h = plot(q(2,:),q(1,:), '.b'); set(h, 'MarkerSize', 25);
    end
    h = plot(x0(2,:),x0(1,:), '.r'); set(h, 'MarkerSize', 25);
    saveas(gcf, [rep 'dual-propagation-' num2str(i) '.eps'], 'epsc');
end
%CMT


