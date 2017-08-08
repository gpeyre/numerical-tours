%% Logistic Classification
% This tour details the logistic classification method (for 2 classes and
% multi-classes).

perform_toolbox_installation('general');

%% Dataset Loading

%%
% Helpers.

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);

%%
% Load the dataset.

name = 'digits';
load(['ml-' name]);

%%
% Randomly permute it.

A = A(randperm(size(A,1)),:);

%%
% Separate the features \(X\) from the data \(y\) to predict information.

X = A(:,1:end-1);
y = A(:,end);
y = y-min(y)+1;
k = max(y);

%%
% \(p\) is the number of samples, \(n\) is the dimensionality of the features,

[p,n] = size(X);

%%
% Display a few samples digits

q = 5;
clf;
for i=1:k
    I = find(y==i);
    for j=1:q
        f = reshape(X(I(j),:), sqrt(n)*[1 1])';
        subplot(q,k, (j-1)*k+i );
        imagesc(-f); axis image; axis off;
    end
end
colormap gray(256);

%% Dimenionality Reduction and PCA
% In order to display in 2D or 3D the data, dimensionality is needed.
% The simplest method is the principal component analysis, which perform an
% orthogonal linear projection on the principal axsis (eigenvector) of the
% covariance matrix. 

%%
% Compute empirical mean \(m \in \RR^n\) and covariance \(C \in \RR^{n \times n}\).

m = mean(X,1);
Xm = X-repmat(m, [p 1]);
C = Xm'*Xm;

%%
% Display the covariance matrix.

clf;
imagesc(C);

%%
% Compute PCA ortho-basis. 

[U,D] = eig(C); 
[d,I] = sort(diag(D), 'descend');
U = U(:,I);

%%
% Compute the feature in the PCA basis.

z = (U'*Xm')';

%%
% Plot sqrt of the eigenvalues.

clf; 
plot(sqrt(max(d,0)), '.-', 'LineWidth', 2, 'MarkerSize', 30);
axis tight;
SetAR(1/2);

%%
% Display in 2D.

col = {'b' 'g' 'r' 'c' 'm' 'y' 'k'};
ms = 25;
clf; hold on;
lgd = {};
for i=1:min(k,length(col))
    I = find(y==i);
    plot(z(I,1), z(I,2), '.', 'Color', col{i}, 'MarkerSize', ms);
    lgd{end+1} = num2str(i);
end
axis tight; axis equal; box on;
legend(lgd);
SetAR(1);

%% Two Classes Logistic Classification

%%
% _Warning:_ Logisitic classification is actually called "logistic
% regression" in the literature, but it is in fact a classification method.

%%
% Select only two classes.

c = [5, 7]; % selected classes
Xsvg = X; ysvg = y; 
X = []; y = [];
for i=1:2
    I = find(ysvg==c(i));
    X = [X; Xsvg(I,:)];
    y = [y; ones(length(I),1)*c(i)];
end

%% Multi-Classes Logistic Classification