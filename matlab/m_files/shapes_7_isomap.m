%% Manifold Learning with Isomap
% This tour explores the Isomap algorithm for manifold learning.

%%
% The <http://waldron.stanford.edu/~isomap/ Isomap> algorithm is introduced in

%%
% _A Global Geometric Framework for Nonlinear Dimensionality Reduction,_
% J. B. Tenenbaum, V. de Silva and J. C. Langford,
% Science 290 (5500): 2319-2323, 22 December 2000.

%CMT
rep = 'results/shapes_isomap/';
if not(exist(rep))
    mkdir(rep);
end
%CMT

perform_toolbox_installation('signal', 'general', 'graph');

%CMT
rep = 'results/fastmarching_bendinginv_2d/';
if not(exist(rep))
    mkdir(rep);
end
%CMT

%% Graph Approximation of Manifolds
% Manifold learning consist in approximating the parameterization of a
% manifold represented as a point cloud.

%%
% First we load a simple 3D point cloud, the famous Swiss Roll.

%% 
% Number of points.

n = 1000;

%%
% Random position on the parameteric domain.

x = rand(2,n);

%%
% Mapping on the manifold.

v = 3*pi/2 * (.1 + 2*x(1,:));
X  = zeros(3,n);
X(2,:) = 20 * x(2,:);
X(1,:) = - cos( v ) .* v;
X(3,:) = sin( v ) .* v;

%%
% Parameter for display.

ms = 50;
lw = 1.5;
v1 = -15; v2 = 20;

%%
% Display the point cloud.

clf;
scatter3(X(1,:),X(2,:),X(3,:),ms,v, 'filled'); 
colormap jet(256);
view(v1,v2); axis('equal'); axis('off');

%CMT
saveas(gcf, [rep 'isomap-cloud.eps'], 'epsc');
%CMT

%%
% Compute the pairwise Euclidean distance matrix.

D1 = repmat(sum(X.^2,1),n,1);
D1 = sqrt(D1 + D1' - 2*X'*X);

%%
% Number of NN for the graph.

k = 6;

%%
% Compute the k-NN connectivity.

[DNN,NN] = sort(D1);
NN = NN(2:k+1,:);
DNN = DNN(2:k+1,:);

%%
% Adjacency matrix, and weighted adjacency.

B = repmat(1:n, [k 1]);
A = sparse(B(:), NN(:), ones(k*n,1));

%%
% Weighted adjacency (the metric on the graph).

W = sparse(B(:), NN(:), DNN(:));

%%
% Display the graph.

options.lw = lw;
options.ps = 0.01;
clf; hold on;
scatter3(X(1,:),X(2,:),X(3,:),ms,v, 'filled'); 
plot_graph(A, X, options);
colormap jet(256);
view(v1,v2); axis('equal'); axis('off');
zoom(.8);

%CMT
saveas(gcf, [rep 'isomap-graph.eps'], 'epsc');
%CMT

%% Floyd Algorithm to Compute Pairwise Geodesic Distances
% A simple algorithm to compute the geodesic distances between all pairs of
% points on a graph is Floyd iterative algorithm. Its complexity is
% |O(n^3)| where |n| is the number of points. It is thus quite slow for
% sparse graph, where Dijkstra runs in |O(log(n)*n^2)|.

%% 
% Floyd algorithm iterates the following update rule, for |k=1,...,n|

%%
% |D(i,j) <- min(D(i,j), D(i,k)+D(k,j)|,

%%
% with the initialization |D(i,j)=W(i,j)| if |W(i,j)>0|, and
% |D(i,j)=Inf| if |W(i,j)=0|.

%%
% Make the graph symmetric.

D = full(W);
D = (D+D')/2;

%% 
% Initialize the matrix.

D(D==0) = Inf;

%%
% Add connexion between a point and itself.

D = D - diag(diag(D));

%EXO
%% Implement the Floyd algorithm to compute the full distance matrix 
%% |D|, where |D(i,j)| is the geodesic distance between 
for i=1:n
    % progressbar(i,n);
    D = min(D,repmat(D(:,i),[1 n])+repmat(D(i,:),[n 1])); 
end
%EXO

%% 
% Find index of vertices that are not connected to the main manifold.

Iremove = find(D(:,1)==Inf);

%%
% Remove Inf remaining values (disconnected comonents).

D(D==Inf) = 0;

%% Isomap with Classical Multidimensional Scaling
% Isomap perform the dimensionality reduction by applying multidimensional
% scaling.

%%
% Please refers to the tours on Bending Invariant for detail on 
% Classical MDS (strain minimization).

%EXO
%% Perform classical MDS to compute the 2D flattening.
% centered kernel
J = eye(n) - ones(n)/n;
K = -1/2 * J*(D.^2)*J;
% diagonalization
opt.disp = 0; 
[Xstrain, val] = eigs(K, 2, 'LR', opt);
Xstrain = Xstrain .* repmat(sqrt(diag(val))', [n 1]);
Xstrain = Xstrain';
% plot graph
clf; hold on;
scatter(Xstrain(1,:),Xstrain(2,:),ms,v, 'filled'); 
plot_graph(A, Xstrain, options);
colormap jet(256);
axis('equal'); axis('off'); 
%EXO

%%
% Redess the points using the two leading eigenvectors of the covariance
% matrix (PCA correction).

[U,L] = eig(Xstrain*Xstrain' / n);
Xstrain1 = U'*Xstrain;

%%
% Remove problematic points.

Xstrain1(:,Iremove) = Inf;

%%
% Display the final result of the dimensionality reduction.


clf; hold on;
scatter(Xstrain1(1,:),Xstrain1(2,:),ms,v, 'filled'); 
plot_graph(A, Xstrain1, options);
colormap jet(256);
axis('equal'); axis('off'); 

%CMT
saveas(gcf, [rep 'isomap-strain.eps'], 'epsc');
%CMT

%%
% For comparison, the ideal locations on the parameter domain.

Y = cat(1, v, X(2,:));
Y(1,:) = rescale(Y(1,:), min(Xstrain(1,:)), max(Xstrain(1,:)));
Y(2,:) = rescale(Y(2,:), min(Xstrain(2,:)), max(Xstrain(2,:)));

%%
% Display the ideal graph on the reduced parameter domain.

clf; hold on;
scatter(Y(1,:),Y(2,:),ms,v, 'filled'); 
plot_graph(A,  Y, options);
colormap jet(256);
axis('equal'); axis('off'); 
camroll(90);

%CMT
saveas(gcf, [rep 'isomap-ideal.eps'], 'epsc');
%CMT


%% Isomap with SMACOF Multidimensional Scaling
% It is possible to use SMACOF instead of classical scaling.

%%
% Please refers to the tours on Bending Invariant for detail on both
% Classical MDS (strain minimization) and SMACOF MDS (stress minimization).

%EXO
%% Perform stress minimization MDS using SMACOF to compute the 2D flattening.
niter = 150;
stress = [];
Xstress = X;
ndisp = [1 5 10 min(niter,100) Inf];
k = 1;
clf;
for i=1:niter
    if ndisp(k)==i
        subplot(2,2,k);
        hold on;
        scatter3(Xstress(1,:),Xstress(2,:),Xstress(3,:),ms,v, 'filled');
        plot_graph(A, Xstress, options);
        colormap jet(256);
        view(v1,v2); axis('equal'); axis('off');
        k = k+1;
    end
    % Compute the distance matrix.
    D1 = repmat(sum(Xstress.^2,1),n,1);
    D1 = sqrt(D1 + D1' - 2*Xstress'*Xstress);
    % Compute the scaling matrix.
    B = -D./max(D1,1e-10);
    B = B - diag(sum(B));
    % update
    Xstress = (B*Xstress')' / n;
    % Xstress = Xstress-repmat(mean(Xstress,2), [1 n]);
    % record stress
    stress(end+1) = sqrt( sum( abs(D(:)-D1(:)).^2 ) / n^2 );
end
%EXO


%%
% Plot stress evolution during minimization.

clf;
plot(stress(1:end), '.-');
axis('tight');

%CMT
%% Perform stress minimization MDS using SMACOF to compute the 2D flattening.
niter = 100;
Xstress = X;
ndisp = [1 5 10 niter Inf];
k = 1;
clf;
for i=1:niter
    if ndisp(k)==i
        clf; hold on;
        scatter3(Xstress(1,:),Xstress(2,:),Xstress(3,:),ms,v, 'filled');
        plot_graph(A, Xstress, options);
        colormap jet(256);
        view(v1,v2); axis('equal'); axis('off');
        zoom(.8);
        saveas(gcf, [rep 'isomap-stress-' num2str(k) '.eps'], 'epsc');
        k = k+1;
    end
    % Compute the distance matrix.
    D1 = repmat(sum(Xstress.^2,1),n,1);
    D1 = sqrt(D1 + D1' - 2*Xstress'*Xstress);
    % Compute the scaling matrix.
    B = -D./max(D1,1e-10);
    B = B - diag(sum(B));
    % update
    Xstress = (B*Xstress')' / n;
end
%CMT

%%
% Compute the main direction of the point clouds.

[U,L] = eig(Xstress*Xstress' / n);
[L,I] = sort(diag(L));
U = U(:,I(2:3));

%%
% Project the points on the two leading eigenvectors of the covariance
% matrix (PCA final projection).

Xstress1 = U'*Xstress;

%%
% Remove problematic points.

Xstress1(:,Iremove) = Inf;

%%
% Display the final result of the dimensionality reduction.

clf; hold on;
scatter(Xstress1(1,:),Xstress1(2,:),ms,v, 'filled'); 
plot_graph(A, Xstress1, options);
colormap jet(256);
axis('equal'); axis('off'); 


%CMT
saveas(gcf, [rep 'isomap-stress.eps'], 'epsc');
%CMT


%% Learning Manifold of Patches
% Isomap algorithm can be used to analyze the structure of a high
% dimensional library of images.

%EXO
%% Apply Isomap to a library of small images, for instance binary digits or
%% faces with a rotating camera.
% No correction for this exercise.
%EXO