%% Linear Regression and Kernel Methods
% This tour studies linear regression method, and its non-linear variant
% using kernlization.

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
% We test the method on the
% <http://scikit-learn.org/stable/modules/generated/sklearn.datasets.load_boston.html
% Boston house prices dataset>, consisting in \(n=506\) samples with
% features \(x_i \in \RR^p\) in dimension \(p=13\). The goal is to predict the price value
% \(y_i \in \RR\). 

%%
% Helpers.

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', 20);
Xm = @(X)X-repmat(mean(X,1), [size(X,1) 1]);
Cov = @(X)Xm(X)'*Xm(X);

%%
% Load the dataset.

name = 'boston_house_prices';
load(['ml-' name]);

%%
% Randomly permute it.

A = A(randperm(size(A,1)),:);

%%
% Separate the features \(X\) from the data \(y\) to predict information.

X = A(:,1:end-1);
y = A(:,end);

%%
% \(n\) is the number of samples, \(p\) is the dimensionality of the features,

[n,p] = size(X);

%% Dimenionality Reduction and PCA
% In order to display in 2-D or 3-D the data, dimensionality is needed.
% The simplest method is the principal component analysis, which perform an
% orthogonal linear projection on the principal axsis (eigenvector) of the
% covariance matrix.

%%
% Display the covariance matrix.

clf;
imagesc(Cov(X));

%%
% Compute PCA ortho-basis and 
% the feature in the PCA basis.

[U,D,V] = svd(Xm(X),'econ');
Z = Xm(X) * V;


%%
% Plot sqrt of the eigenvalues.

clf; 
plot(diag(D), '.-', 'LineWidth', 2, 'MarkerSize', 30);
axis tight;
SetAR(1/2);

%%
% Display the points cloud of feature vectors in 3-D PCA space.

options.disp_dim = 3;
clf; plot_multiclasses(X,ones(n,1),options);
SetAR(1);

%%
% 1D plot of the function to regress along the main eigenvector axes.

col = {'b' 'g' 'r' 'c' 'm' 'y' 'k'};
clf;
for i=1:min(p,3)
    subplot(3,1,i);
    plot(Z(:,i), y, '.', 'Color', col{i}, 'MarkerSize', 20);
    axis tight;
end

%% Linear Regression
% We look for a linear relationship 
%   \( y_i = \dotp{w}{x_i} \)
% written in matrix format 
%   \( y= X w \)
% where the rows of \(X \in \RR^{n \times p}\) stores the features \(x_i \in \RR^p\).

%%
% Since here \( n > p \), this is an over-determined system, which can
% solved in the least square sense
%   \[ \umin{ w }  \norm{Xw-y}^2 \]
% whose solution is given using the Moore-Penrose pseudo-inverse
%   \[ w = (X^\top X)^{-1} X^\top y \]


%%
% Split into training and testing.

n0 = round(.5*n); n1 = n-n0;
X0 = X(1:n0,:);     y0 = y(1:n0);
X1 = X(n0+1:end,:); y1 = y(n0+1:end);

%%
% Least square solution.

w = (X0'*X0) \ (X0'*y0);

%%
% Mean-square error on testing set.

E = sqrt( sum( (X1*w-y1).^2 ) / n1 );

%%
% Regularization is obtained by introducing a penalty. It is often called
% ridge regression, and is defined as 
%   \[ \umin{ w }  \norm{Xw-y}^2 + \lambda \norm{w}^2 \]
% where \(\lambda>0\) is the regularization parameter.

%%
% The solution is given using the following equivalent formula
%   \[ w = (X^\top X + \lambda \text{Id}_p )^{-1} X^\top y, \]
%   \[ w = X^\top ( XX^\top + \lambda \text{Id}_n)^{-1} y, \]
% When \(p<n\) (which is the case here), the first formula should be
% prefered. 

%%
% In contrast, when the dimensionality \(p\) of the feature is very
% large and there is little data, the second is faster. Furthermore, this
% second expression is generalizable to Kernel Hilbert space setting,
% corresponding possibly to \(p=+\infty\) for some kernels. 

lambda = .1;
w = (X0'*X0+lambda*eye(p)) \ (X0'*y0);
w1 = X0'*( (X0*X0'+lambda*eye(n0)) \ y0 );
fprintf('Error (should be 0): %.4f\n', norm(w-w1)/norm(w));

%EXO
%% Display the evolution of the error \(E\) as a function of \(\lambda\).
q = 50;
lambda_list = linspace(0,20,q);
W = []; E = [];
for i=1:q
    lambda = lambda_list(i);
    w = (X0'*X0+lambda*eye(p)) \ (X0'*y0);
    W(:,i) = w; % bookkeeping
    E(i) = sqrt( sum( (X1*w-y1).^2 ) / n1 );
end
% Display error evolution.
clf;
plot(lambda_list, E, 'LineWidth', 2);
axis tight;
SetAR(1/2);
xlabel('\lambda');
ylabel('E');
%EXO

%EXO
%% Display the regularization path, i.e. the evolution of \(w\) as a function 
%% of \(\lambda\).
clf;
plot(lambda_list, W', 'LineWidth', 2);
axis tight;
SetAR(1/2);
xlabel('\lambda'); ylabel('w_i');
%EXO

%EXO
%% Perform feature selection using \(\ell^1\) regularization (aka the Lasso), 
%%   \[ \umin{ w }  \frac{1}{2}\norm{Xw-y}^2 + \lambda \norm{w}_1. \]
%EXO


%% Kernelized Ridge Regression
% In order to perform non-linear and non-parametric regression, it is
% possible to use  kernelization. It is non-parametric in the sense that the number of
% parameter grows with the number \(n\) of samples (while for the initial
% linear  method, the number of parameter is \(p\)). This allows in particular to
% generate estimator of arbitrary complexity. 

%%
% Given a kernel \( \kappa(x,z) \in \RR \) defined for \((x,z) \in \RR^p \times \RR^p\),
% the kernelized method replace the linear approximation functional \(f(x) =
% \dotp{x}{w}\) by a sum of kernel centered on the samples 
% \[ f_h(x) = \sum_{i=1}^n h_i k(x_i,x) \]
% where \(h \in \RR^n\) is the unknown vector of weight to find. 

%%
% When using the linear kernel \(\kappa(x,y)=\dotp{x}{y}\), one retrieves
% the previously studied linear method. 

%%
% Generate synthetic data in 2D.
% Add noise to a deterministic map.

B = 3;
n = 500; p = 2;
X = 2*B*rand(n,2)-B;
rho = .5; % noise level
y = peaks(X(:,1), X(:,2)) + randn(n,1)*rho;

%%
% Display as scattered plot.

clf;
scatter(X(:,1), X(:,2), ones(n,1)*20, y, 'filled');
colormap jet(256); 
axis equal; axis([-B B -B B]); box on;

%% 
% Macro to compute pairwise squared Euclidean distance matrix.

distmat = @(X,Z)bsxfun(@plus,dot(X',X',1)',dot(Z',Z',1))-2*(X*Z');


%%
% The gaussian kernel is the most well known and used kernel
% \[ \kappa(x,y) \eqdef e^{-\frac{\norm{x-y}^2}{2\sigma^2}} . \]
% The bandwidth parameter \(\si>0\) is crucial and controls the locality of
% the model. It is typically tuned through cross validation.  

sigma = .3;
kappa = @(X,Z)exp( -distmat(X,Z)/(2*sigma^2) );


%%
% Once avaluated on grid points, the kernel define a matrix
% \[ K = (\kappa(x_i,x_j))_{i,j=1}^n \in \RR^{n \times n}.  \]

K = kappa(X,X);

%%
% The weights \(h \in \RR^n \) are solutions of
%   \[ \umin{h} \norm{Kh-y}^2 + \la \dotp{Kh}{h}  \]
% and hence can be computed by solving a linear system
%   \[ h = (K+\la \text{Id}_n)^{-1} y  \]

lambda = 0.01;
h = (K+lambda*eye(n))\y;

%%
% Regressor.

Y = @(x)kappa(x,X)*h;

%%
% Evaluation on a 2D grid.

q = 101;
t = linspace(-B,B,q);
[v,u] = meshgrid(t,t);
Xn = [u(:), v(:)];

%%
% Display as an image.

yn = reshape(Y(Xn),[q,q]);
clf;
imagesc(t,t,yn); axis image; axis off; 
colormap jet(256);
    
%EXO
%% Display the evolution of the regression as a function of \(\sigma\).
sigma_list = [.05 .1 .5 1 5 10];
%
clf;
for i=1:length(sigma_list)
    sigma = sigma_list(i);
    kappa = @(X,Z)exp( -distmat(X,Z)/(2*sigma^2) );
    % Regressor.
    h = (kappa(X,X)+lambda*eye(n))\y;
    Y = @(x)kappa(x,X)*h;
    % Eval on the grid
    yn = reshape(Y(Xn),[q,q]);
    %
    subplot(2,3,i);
    imagesc(t,t,yn); axis image; axis off; 
    colormap jet(256);
    title(['\sigma=' num2str(sigma)]);
end
%EXO

%EXO
%% Apply the kernelize regression to a real life dataset. Study the influence of \(\la\) and \(\si\).
%EXO
