%% Entropic Regularization of Optimal Transport
% This numerical tours exposes the general methodology of regularizing the
% optimal transport (OT) linear program using entropy. This allows to
% derive fast computation algorithm based on iterative projections
% according to a Kulback-Leiber divergence.
% \[ \DeclareMathOperator{\KL}{KL}
% \newcommand{\KLdiv}[2]{\KL\pa{#1 | #2}}
% \newcommand{\KLproj}{P^{\tiny\KL}}
% \def\ones{\mathbb{I}} \]

perform_toolbox_installation('signal', 'general');

%% Entropic Regularization of Optimal Transport
% We consider two input histograms \(p,q \in \Si_N\), where we denote the simplex in \(\RR^N\)
% \[ \Si_{N} \eqdef \enscond{ p \in (\RR^+)^N }{ \sum_i p_i = 1 }.  \]
% We consider the following discrete regularized transport
% \[  W_\ga(p,q) \eqdef \umin{\pi \in \Pi(p,q)} \dotp{C}{\pi} - \ga E(\pi).  \]
% where the polytope of coupling is defined as
% \[ \Pi(p,q) \eqdef \enscond{\pi \in (\RR^+)^{N \times N}}{ \pi \ones = p, \pi^T \ones = q },  \]
% where \( \ones \eqdef (1,\ldots,1)^T \in \RR^N \),
% and for \(\pi \in (\RR^+)^{N \times N}\), we define its entropy as
% \[ E(\pi) \eqdef - \sum_{i,j} \pi_{i,j} ( \log(\pi_{i,j}) - 1). \]

%%
% When \(\ga=0\) one recovers the classical (discrete) optimal transport.
% We refer to the monograph <#biblio [Villani]> for more details about OT.
% The idea of regularizing transport to allows for faster computation is
% introduced in <#biblio [Cuturi]>.

%%
% Here the matrix \(C \in (\RR^+)^{N \times N} \) defines the ground cost, i.e.
% \(C_{i,j}\) is the cost of moving mass from a bin indexed by \(i\) to a bin indexed by \(j\).

%%
% The regularized transportation problem can be re-written as a projection
% \[ W_\ga(p,q) = \ga�\umin{\pi \in \Pi(p,q)} \KLdiv{\pi}{\xi}
% 	\qwhereq
% 	\xi_{i,j} = e^{ -\frac{C_{i,j}}{\ga} }  \]
% of \(\xi\) according to the Kullback-Leibler divergence.
% The Kullback-Leibler divergence between \(\pi, \xi \in (\RR^+)^P\) is
% \[ \KLdiv{\pi}{\xi} = \sum_{i,j} \pi_{i,j} \pa{�\log\pa{ \frac{\pi_{i,j}}{\xi_{i,j}} } - 1}. \]

%%
% This interpretation of regularized transport as a KL projection and its numerical
% applications are detailed in <#biblio [BenamouEtAl]>.

%%
% Given a convex set \(\Cc \subset \RR^N\), the projection according to the Kullback-Leiber divergence is defined as
% \[ \KLproj_\Cc(\xi) = \uargmin{ \pi \in \Cc } \KLdiv{\pi}{\xi}. \]

%% Iterative Bregman Projection Algorithm
% Given affine constraint sets \( (\Cc_1,\ldots,\Cc_K) \), we aim at computing
% \[   \KLproj_\Cc(\xi) \qwhereq \Cc = \Cc_1 \cap \ldots \cap \Cc_K. \]

%%
% This can be achieved, starting by \(\pi_0=\xi\), by iterating
% \[ \forall \ell \geq 0, \quad \pi_{\ell+1} =  \KLproj_{\Cc_\ell}(\pi_\ell), \]
% where the index of the constraints should be understood modulo \(K\),
% i.e. we set \( \Cc_{\ell+K}=\Cc_\ell \).

%%
% One can indeed show that \(\pi_\ell \rightarrow \KLproj_\Cc(\bar \pi)\).
% We refer to <#biblio [BauschkeLewis]> for more details about this
% algorithm and its extension to compute the projection on the intersection of
% convex sets (Dikstra algorithm).

%% Iterative Projection for Regularized Transport aka Sinkhorn's Algorithm
% We can re-cast the regularized optimal transport problem within this
% framework by introducing
% \[ \Cc_1 \eqdef \enscond{\pi \in (\RR^+)^{N \times N} }{\pi \ones = p}
% \qandq
%  \Cc_2 \eqdef \enscond{\pi \in (\RR^+)^{N \times N} }{\pi^T \ones = q}. \]

%%
% The KL projection on \(\Cc_1\) sets are easily computed by divisive
% normalization of rows. Indeed, denoting
% \( \pi = \KLproj_{\Cc_1}(\bar \pi) \), one has
% \[ \forall (i,j), \quad
%   \pi_{i,j} = \frac{ p_i \bar\pi_{i,j} }{ \sum_{s} \bar\pi_{i,s} } \]
% and similarely for \(\KLproj_{\Cc_1}(\bar \pi) \) by replacing rows by
% colums.

%%
% A fundamental remark is that, if \(\bar\pi = \diag(a) \xi \diag(b)\) (a
% so-called diagonal scaling of the kernel \(\xi\)), then one has
% \[ \KLproj_{\Cc_1}(\bar \pi) = \diag(\tilde a) \xi \diag(b)
%   \qandq
%   \KLproj_{\Cc_2}(\bar \pi) = \diag(a) \xi \diag(\tilde b)\]
%  where the new scaling reads
%   \[ \tilde a = \frac{p}{\xi(b)}  \qandq  \tilde b = \frac{q}{\xi^T(a)} \]
% where \(\frac{\cdot}{\cdot}\) is entry-wise division.

%%
% This means that the iterates of Bregman iterative projections, starting
% with \( a_0 \eqdef \ones \) always have the form
% \( \pi_\ell = \diag(a_\ell) \xi \diag(b_\ell)  \)
% and these diagonal scaling weights are updated as follow
% \[  a_{\ell+1} \eqdef \frac{p}{\xi(b_\ell)}
%       \qandq
%   b_{\ell+1} \eqdef \frac{q}{\xi^T(a_{\ell+1})}. \]
% This algorithm is in fact the well known Sinkhorn algorithm <#biblio [Sinkhorn]>.

%% Transport Between Point Clouds
% We first test the method for two input measures that are uniform measures
% (i.e. constant histograms) supported on two point clouds
% (that do not necessarily have the same size).

%%
% We thus first load two points clouds \(x=(x_i)_{i=1}^{N_1}, y=(y_i)_{i=1}^{N_2}, \)
% where \(x_i, y_i \in \RR^2\).

%%
% Number of points in each cloud.

N = [300,200];

%%
% Dimension of the clouds.

d = 2;

%%
% Point cloud \(x\), of \(N_1\) points inside a square.

x = rand(2,N(1))-.5;

%%
% Point cloud \(y\), of \(N_2\) points inside an anulus.

theta = 2*pi*rand(1,N(2));
r = .8 + .2*rand(1,N(2));
y = [cos(theta).*r; sin(theta).*r];

%%
% Shortcut for displaying point clouds.

plotp = @(x,col)plot(x(1,:)', x(2,:)', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 2);

%%
% Display of the two clouds.

clf; hold on;
plotp(x, 'b');
plotp(y, 'r');
axis('off'); axis('equal');

%%
% Cost matrix \(C_{i,j} = \norm{x_i-y_j}^2\).

x2 = sum(x.^2,1); y2 = sum(y.^2,1);
C = repmat(y2,N(1),1)+repmat(x2.',1,N(2))-2*x.'*y;

%%
% Target histograms, here uniform histograms.

p = ones(N(1),1)/N(1);
q = ones(N(2),1)/N(2);

%%
% Regularization strength \(\ga\).

gamma = .01;

%%
% Gibbs Kernel.

xi = exp(-C/gamma);

%%
% Initialization of \(b_0=\ones_{N_2}\) (\(a_0\) does not need to be
% initialized).

b = ones(N(2),1);

%%
% One sinkhorn iterations.

a = p ./ (xi*b);
b = q ./ (xi'*a);

%EXO
%% Implement Sinkhorn algorithm.
%% Display the evolution of the constraints satisfaction errors
%% \( \norm{ \pi \ones - p }, \norm{ \pi^T \ones - q } \)
%% (you need to think about how to compute these residuals from \((a,b)\) alone).
b = ones(N(2),1);
niter = 300;
Err_p = []; Err_q = [];
for i=1:niter
    a = p ./ (xi*b);
    Err_q(end+1) = norm( b .* (xi'*a) - q )/norm(q);
    b = q ./ (xi'*a);
    Err_p(end+1) = norm( a .* (xi*b) - p )/norm(p);
end
% Display the violation of constraint error in log-plot.
clf;
subplot(2,1,1);
plot(log10(Err_p)); axis tight; title('log|\pi 1 - p|');
subplot(2,1,2);
plot(log10(Err_q)); axis tight; title('log|\pi^T 1 - q|');
%EXO

%%
% Compute the final matrix.

Pi = diag(a)*xi*diag(b);

%%
% Display it.

clf;
imageplot(Pi);

%EXO
%% Display the regularized transport solution for various values of \(\gamma\).
%% For a too small value of \(\gamma\), what do you observe ?
glist = [.1 .01 .005 .001 ];
niter = 300;
clf;
for k=1:length(glist)
    gamma = glist(k);
    xi = exp(-C/gamma);
    b = ones(N(2),1);
    for i=1:niter
        a = p ./ (xi*b);
        b = q ./ (xi'*a);
    end
    Pi = diag(a)*xi*diag(b);
    imageplot( clamp(Pi,0,min(1./N)*.3) , ['\gamma=' num2str(gamma)], 2,2,k);
end
%EXO

%%
% Compute the obtained optimal \(\pi\).

Pi = diag(a)*xi*diag(b);

%%
% Keep only the highest entries of the coupling matrix, and use them to
% draw a map between the two clouds.
% First we draw "strong" connexions, i.e. linkds \((i,j)\) corresponding to
% large values of \(\pi_{i,j}\).
% We then draw weaker connexions.

clf;
hold on;
A = sparse( Pi .* (Pi> min(1./N)*.7) ); [i,j,~] = find(A);
h = plot([x(1,i);y(1,j)], [x(2,i);y(2,j)], 'k');
set(h, 'LineWidth', 2); % weaker connections.
A = sparse( Pi .* (Pi> min(1./N)*.3) ); [i,j,~] = find(A);
h = plot([x(1,i);y(1,j)], [x(2,i);y(2,j)], 'k:');
set(h, 'LineWidth', 1);
plotp(x, 'b'); % plot the two point clouds.
plotp(y, 'r');
axis('off'); axis('equal');

%% Transport Between Histograms
% We now consider a different setup, where the histogram values
% \(p,q\) are not uniform, but the measures are defined on a uniform grid
% \(x_i=y_i=i/N\). They are thue often refered to as "histograms".

%%
% Size \(N\) of the histograms.

N = 200;

%%
% We use here a 1-D square Euclidean metric.

t = (0:N-1)'/N;

%%
% Define the histogram \(p,q\) as translated Gaussians.

Gaussian = @(t0,sigma)exp( -(t-t0).^2/(2*sigma^2) );
normalize = @(p)p/sum(p(:));
sigma = .06;
p = Gaussian(.25,sigma);
q = Gaussian(.8,sigma);

%%
% Add some minimal mass and normalize.

vmin = .02;
p = normalize( p+max(p)*vmin);
q = normalize( q+max(q)*vmin);

%%
% Display the histograms.

clf;
subplot(2,1,1);
bar(t, p, 'k'); axis tight;
subplot(2,1,2);
bar(t, q, 'k'); axis tight;

%%
% Regularization strength \(\ga\).

gamma = (.03)^2;

%%
% The Gibbs kernel is a Gaussian convolution,
% \[ \xi_{i,j} = e^{ -(i/N-j/N)^2/\gamma }. \]

[Y,X] = meshgrid(t,t);
xi = exp( -(X-Y).^2 / gamma);

%%
% Initialization of \(b_0=\ones_{N}\).

b = ones(N,1);

%%
% One sinkhorn iteration.

a = p ./ (xi*b);
b = q ./ (xi'*a);

%EXO
%% Implement Sinkhorn algorithm.
%% Display the evolution of the constraints satisfaction errors
%% \( \norm{ \pi \ones - p }, \norm{ \pi^T \ones - q } \)
%% (you need to think how to compute it from \((a,b)\).
b = ones(N,1);
niter = 2000;
Err_p = []; Err_q = [];
for i=1:niter
    a = p ./ (xi*b);
    Err_q(end+1) = norm( b .* (xi*a) - q )/norm(q);
    b = q ./ (xi'*a);
    Err_p(end+1) = norm( a .* (xi'*b) - p )/norm(p);
end
% Display the violation of constraint error in log-plot.
clf;
subplot(2,1,1);
plot(log10(Err_p)); axis tight; title('log|\pi 1 - p|');
subplot(2,1,2);
plot(log10(Err_q)); axis tight; title('log|\pi^T 1 - q|');
%EXO

%%
% Display the coupling. Use a log domain plot to better vizualize it.

Pi = diag(a)*xi*diag(b);
clf;
imageplot(log(Pi+1e-5));

%%
% One can compute an approximation of the transport plan
% between the two measure by computing the so-called barycentric projection map
% \[ t_i \in [0,1] \longmapsto s_j \eqdef \frac{\sum_{j} \pi_{i,j} t_j }{ \sum_{j} \pi_{i,j} }
%    =  \frac{ [a \odot \xi(b \odot t)]_j }{ p_i }. \]
% where \(\odot\) and \(\frac{\cdot}{\cdot}\) are the enry-wise multiplication and division.

%%
% This computation can thus be done using only multiplication with the
% kernel \(\xi\).

s = (xi*(b.*t)) .* a ./ p;

%%
% Display the transport map, super-imposed over the coupling.

clf; hold on;
imagesc(t,t,log(Pi+1e-5)); colormap gray(256);
plot(s,t, 'r', 'LineWidth', 3);
axis image; axis off; axis ij;

%% Wasserstein Barycenters
% Instead of computing transport, we now turn to the problem of computing
% barycenter of measures. A barycenter \(q\) solves
% \[ \umin{q} \sum_{k=1}^K W_\ga(p_k,q)  \]
% where \(\la_k\) are positive weights with \(\sum_k \la_k=1\). This
% follows the definition of barycenters proposed in
% <#biblio [AguehCarlier]>.

%%
% With of the histograms.

N = 70;

%%
% Load input histograms \( (p_k)_{k=1}^K \).

names = {'disk' 'twodisks' 'letter-x' 'letter-z'};
r = [.35 .19 .12*N .12*N];
vmin = .05;
P = [];
for i=1:length(names)
    options.radius = r(i);
    p = load_image(names{i},N, options);
    p = normalize( rescale(p)+vmin );
    P(:,:,i) = p;
end
K = size(P,3);

%%
% Display the input histograms.

a = mat2cell(P, N,N,ones(K,1));
clf;
imageplot(a, '', 2,2);

%%
% In this specific case, the kernel \(\xi\) associated with the
% squared Euclidean norm is a convolution with a Gaussian filter
% \[ \xi_{i,j} = e^{ -\norm{i/N-j/N}^2/\gamma } \]
% where here \((i,j)\) are 2-D indexes.

%%
% The multiplication against the kernel \(\xi(a_\ell)\)
% can now be computed efficiently, using fast convolution methods.
% This
% crucial points was exploited and generalized in <#biblio [SolomonEtAl]>
% to design fast optimal transport algorithm.

%%
% Regularization strength \(\ga\).

gamma = (.04)^2;

%%
% Define the \(\xi\) kernel.
% We use here the fact that the convolution is separable to implement it
% using only 1-D convolution, which further speeds up computations.

n = 41; % width of the convolution kernel
t = linspace(-n/(2*N),n/(2*N),n)';
g = exp(-t.^2 / gamma); g2 = g*g';
xi = @(x)conv2(conv2(x, g, 'same')', g, 'same')';

%%
% Display the application of the \(\xi\) kernel on one of the input histogram.

clf;
imageplot({P(:,:,1) xi(P(:,:,1))});

%%
% Weights for isobarycenter.

lambda = ones(K,1)/K;

%%
% It is shown in <#biblio [BenamouEtAl]> that the problem of Barycenter computation
% boilds down to
% optimizing over couplings \((\pi_k)_k\), and that this can be achieved
% using iterative Bregman projection that defines iterates \((\pi_{k,\ell})_k\).
% These iterates can be written using diagonal scalings
% \( \pi_{k,\ell} \eqdef \diag(a_{k,\ell}) \xi \diag(b_{k,\ell}). \)

%%
% Initialize the scaling factors.

b = ones(N,N,K); a = b;

%%
% The first step of the Bregman projection method corresponds to the
% projection on the fixed marginals constraints \(\pi_k \ones=p_k\). This
% is achieved by updating
% \[ \forall k=1,\ldots,K, \quad a_{k,\ell+1} = \frac{p_k}{ \xi( b_{k,\ell} ) }. \]

for k=1:K
    a(:,:,k) = P(:,:,k) ./ xi(b(:,:,k));
end

%%
% The second step of the Bregman projection method corresponds to the
% projection on the equal marginals constraints \(\forall k, \pi_k^T \ones=q\). This
% is achieved by first computing the target barycenter using a geometric means
% \[ \log(q_\ell) \eqdef \sum_k \lambda_k \log( b_{k,\ell} \odot \xi( a_{k,\ell} ) ). \]


q = zeros(N);
for k=1:K
    q = q + lambda(k) * log( max(1e-19, b(:,:,k) .* xi(a(:,:,k)) ) );
end
q = exp(q);

%%
% And then one can update the other maginals to be equal to this barycenter at step \(\ell\).
% \[ \forall k=1,\ldots,K, \quad b_{k,\ell+1} \eqdef \frac{q_\ell}{ \xi(a_{k,\ell+1}) }. \]

for k=1:K
    b(:,:,k) = q ./ xi(a(:,:,k));
end


%EXO
%% Implement the iterative algorithm to compute the iso-barycenter of the measures.
%% Plot the decay of the error \( \sum_k \norm{\pi_k \ones - p_k} \).
niter = 800;
b = ones(N,N,K); a = b;
mynorm = @(x)norm(x);
Err_q = [];
for i=1:niter
    Err_q(i) = 0;
    for k=1:K
        Err_q(i) = Err_q(i) + mynorm( a(:,:,k) .* xi(b(:,:,k)) - P(:,:,k) )/mynorm(P(:,:,k));
        a(:,:,k) = P(:,:,k) ./ xi(b(:,:,k));
    end
    q = zeros(N);
    for k=1:K
        q = q + lambda(k) * log( max(1e-19, b(:,:,k) .* xi(a(:,:,k)) ) );
    end
    q = exp(q);
    % clf; imageplot(q); drawnow;
    for k=1:K
        b(:,:,k) = q ./ xi(a(:,:,k));
    end
end
clf;
plot(log(Err_q)); axis tight;
%EXO

%%
% Display the barycenter.

clf;
imageplot(q);

%EXO
%% Compute barycenters for varying weights \(\la\) corresponding to
%% a bilinear interpolation inside a square.
m = 5;
[T,S] = meshgrid(linspace(0,1,m), linspace(0,1,m));
clf;
niter = 1000;
for j=1:m^2
    % weights
    lambda = [S(j)*T(j), (1-S(j))*T(j), S(j)*(1-T(j)), (1-S(j))*(1-T(j))];
    % computation
    b = ones(N,N,K); a = b;
    for i=1:niter
        for k=1:K
            a(:,:,k) = P(:,:,k) ./ xi(b(:,:,k));
        end
        q = zeros(N);
        for k=1:K
            q = q + lambda(k) * log( max(1e-19, b(:,:,k) .* xi(a(:,:,k)) ) );
        end
        q = exp(q);
        for k=1:K
            b(:,:,k) = q ./ xi(a(:,:,k));
        end
    end
    % display
    subplot(m,m,j);
    imageplot(q); drawnow;
end
%EXO



%% Bibliography
% <html><a name="biblio"></a></html>

%%
% * [Villani]�C. Villani, (2009). Optimal transport: old and new, volume 338. Springer Verlag.
% * [Cuturi]�M. Cuturi, (2013). Sinkhorn distances: Lightspeed computation of optimal transport. In Burges, C. J. C., Bottou, L., Ghahramani, Z., and Weinberger, K. Q., editors, Proc. NIPS, pages 2292-2300.
% * [AguehCarlier]�M. Agueh, and G Carlier, (2011). Barycenters in the Wasserstein space. SIAM J. on Mathematical Analysis, 43(2):904-924.
% * [CuturiDoucet]�M. Cuturi and A. Doucet (2014). Fast computation of wasserstein barycenters. In Proc. ICML.
% * [BauschkeLewis]�H. H. Bauschke and A. S. Lewis. Dykstra's algorithm with Bregman projections: a convergence proof. Optimization, 48(4):409-427, 2000.
% * [Sinkhorn] R. Sinkhorn. A relationship between arbitrary positive matrices and doubly stochastic matrices. Ann. Math. Statist., 35:876-879, 1964.
% * [SolomonEtAl] J. Solomon, F. de Goes, G. Peyr�, M. Cuturi, A. Butscher, A. Nguyen, T. Du, and L. Guibas. Convolutional Wasserstein distances: Efficient optimal transportation on geometric domains. Transaction on Graphics, Proc. SIGGRAPH, 2015.
% * [BenamouEtAl] J-D. Benamou, G. Carlier, M. Cuturi, L. Nenna, G. Peyr�. Iterative Bregman Projections for Regularized Transportation Problems. SIAM Journal on Scientific Computing, 37(2), pp. A1111-A1138, 2015.
