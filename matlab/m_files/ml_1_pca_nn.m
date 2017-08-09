%% PCA, Nearest-Neighbors and Clustering
% This tour details PCA (visualization), NN classification and K-means on the IRIS dataset.

%%
% We recommend that after doing this Numerical Tours, you apply it to your
% own data, for instance using a dataset from <https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/ LibSVM>.

%%
% _Disclaimer:_ these machine learning tours are intended to be
% overly-simplistic implementations and applications of baseline machine learning methods. 
% For more advanced uses and implementations, we recommend
% to use a state-of-the-art library, the most well known being
% <http://scikit-learn.org/ Scikit-Learn>

perform_toolbox_installation('general');

%% Dataset Loading

%%
% Helpers.

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

%%
% Load the dataset.

name = 'iris';
load(['ml-' name]);

%%
% Randomly permute it.

A = A(randperm(size(A,1)),:);


%%
% Separate the features from the class information.
% Be sure to start the class at index 1. 

X = A(:,1:end-1);
y = A(:,end);
y = y-min(y)+1;

%%
% \(p\) is the number of samples, \(n\) is the dimensionality of the features,
% \(k\) is the number of classes.

[p,n] = size(X);
k = max(y);

%% Dimenionality Reduction and PCA
% In order to display in 2D or 3D the data, dimensionality is needed.
% The simplest method is the principal component analysis, which perform an
% orthogonal linear projection on the principal axsis (eigenvector) of the
% covariance matrix.

%%
% Compute empirical mean \(m \in \RR^n\) and covariance \(C \in \RR^{n \times n}\).

Xm = @(X)X-repmat(mean(X,1), [size(X,1) 1]);
Cov = @(X)Xm(X)'*Xm(X);

%%
% Display the covariance matrix.

clf;
imagesc(Cov(X));

%%
% Compute PCA ortho-basis. 

[U,D,V] = svd(Xm(X),'econ');

%%
% Compute the feature in the PCA basis.

Z = Xm(X) * V;

%%
% Plot sqrt of the eigenvalues.

clf; 
plot(diag(D), '.-', 'LineWidth', 2, 'MarkerSize', 30);
axis tight;
SetAR(1/2);


%%
% Display in 2D.

col = {'b' 'g' 'r' 'c' 'm' 'y' 'k'};
ms = 25;
clf; hold on;
for i=1:k
    I = find(y==i);
    plot(Z(I,1), Z(I,2), '.', 'Color', col{i}, 'MarkerSize', ms);
end
axis tight; axis equal; box on;
SetAR(1);


%%
% Display in 2D.

clf; hold on;
for i=1:k
    I = find(y==i);
    plot3(Z(I,1), Z(I,2), Z(I,3), '.', 'Color', col{i}, 'MarkerSize', ms);
end
view(3); axis tight; axis equal; box on;
SetAR(1);

%% Supervised Learning: Nearest Neighbor Classification

%%
% Split into training and testing.

p0 = round(.5*p);
p1 = p-p0;

X0 = X(1:p0,:);     y0 = y(1:p0);
X1 = X(p0+1:end,:); y1 = y(p0+1:end);

%% 
% Macro to compute pairwise squared Euclidean distance matrix.

distmat = @(X,Z)bsxfun(@plus,dot(X',X',1)',dot(Z',Z',1))-2*(X*Z');

%%
% Compute Euclidean distance between some \(x_{1,i}\) (for a fixed \(i\)) in the testing set
% and all other \(x_{1,j}\) in the training set. 

i = 1;
D = distmat(X0,X1(i,:));

%%
% Sort the distance and generate the list of sorted classes \(Y\).

[~,I] = sort(D);
Y = y(I);

%%
% Display class evolution as distance grows.

clf;
plot(Y, '.', 'MarkerSize', ms); 
axis tight; box on;
SetAR(1);

%%
% Perform a \(R\)-neareast neighbor classification.
% The output class \(C\) is the one which is the most represented among
% the \(R\) closest points. 

R = 5;  % number of NN
h = hist(Y(1:R,:), 1:k);
[~,C] = max(h);

%EXO
%% Do the same, but for all the point in the test set, and for varying \(R\).
%% Show how the classification score \(S\) (number of correctly classified points)
%% evolves with \(R\)
D = distmat(X0,X1);
[Ds,I] = sort(D,1);
Y = y(I);
Rmax = 50;
for R=1:Rmax
    if R==1
        C = Y(1,:); 
    else
        h = hist(Y(1:R,:), 1:k);
        [~,C] = max(h);
    end
    % correct classification 
    S(R) = sum(C(:)==y1)/p1;
end
clf;
plot(1:Rmax, S, '.-', 'MarkerSize', ms);
axis tight;
SetAR(1/2);
xlabel('R'); ylabel('S');
%EXO

%EXO
%% Display, as a function of the position in 2D PCA space, the class output by 
%% the R-NN method.
% bounding boxes
B = max(max(abs(Z(:,1:2))));
q = 200;
r = linspace(-B,B,q);
[V,U] = meshgrid(r,r);
z1 = [U(:),V(:)];
% test for different R
Rlist = [1 5 10 40];
clf;
for i=1:length(Rlist)
    R=Rlist(i);
    %
    D = distmat(Z(:,1:2),z1);
    [Ds,I] = sort(D,1);
    Y = y(I);
    %
    if R==1
        C = Y(1,:);
    else
        h = hist(Y(1:R,:), 1:k);
        [~,C] = max(h);
    end
    C = reshape(C, [q q]);
    % display
    subplot(2,2,i);
    hold on;
    imagesc(r,r,C');
    for i=1:k
        I = find(y==i);
        plot(Z(I,1), Z(I,2), '.', 'Color', col{i}, 'MarkerSize', ms);
    end
    axis tight; axis equal; axis off;
    title(['R=' num2str(R)]);
    SetAR(1);
end
colormap jet(256);
%EXO

%% Unsupervised Learning: K-means



