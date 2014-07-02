%% Subdivision Curves
% Subdvision methods progressively refine a discrete curve and
% converge to a smooth curve. This allows to perform an
% interpolation or approximation of a given coarse dataset.

perform_toolbox_installation('signal', 'general', 'graph', 'wavelet_meshes');

%% Curve Subdivision 
% Starting from an initial set of control points (which can be seen as a
% coarse curve), subdivision produces a smooth 2-D curve.

%%
% Shortcut to plot a periodic curve.

ms = 20; lw = 1.5;
myplot = @(f,c)plot(f([1:end 1]), c, 'LineWidth', lw, 'MarkerSize', ms);
myaxis = @(rho)axis([-rho 1+rho -rho 1+rho], 'off');

%%
% We represent a dicretized curve of \(N\) points as a vector of complex numbers \(f \in \CC^N\).
% Since we consider periodic boundary conditions, we assume the vectors
% have periodic boundary conditions. 

%%
% Define the initial coarse set of control points \(f_0 \in \CC^{N_0}\).

f0 =    [0.11 0.18 0.26 0.36 0.59 0.64 0.80 0.89 0.58 0.22 0.18 0.30 0.58 0.43 0.42]' + ...
   1i * [0.91 0.55 0.91 0.58 0.78 0.51 0.81 0.56 0.10 0.16 0.35 0.42 0.40 0.24 0.31]';
  
%%
% Rescale it to fit in a box.

f0 = rescale(real(f0),.01,.99) + 1i * rescale(imag(f0),.01,.99);

%%
% Display it.

clf; myplot(f0, 'k.-'); 
myaxis(0);


%%
% One subdivision step reads
% \[ f_{j+1} = (f_j \uparrow 2) \star h. \]
% This produces discrete curves \(f_j \in \CC^{N_j}\) where \(N_j = N_0
% 2^j\). 

%%
% Here \(\uparrow 2\) is the up-sampling operator
% \[ (f \uparrow 2)_{2i}=f_i \qandq (f \uparrow 2)_{2i+1} = 0.  \]

%%
% Recall that the periodic discrete convolution is defined as
% \[ (f \star h)_i = \sum_j f_j h_{i-j}, \]
% where the filter \(h\) is zero-padded to reach the same size as \(f\).

%%
% The low pass filter (subdivision kernel) \(h \in \CC^K\) should
% satisfies 
% \[ \sum_i h_i = 2 . \]
% This ensure that the center of gravity of the curve stays constant
% \[ \frac{1}{N_j} \sum_{i=1}^{N_j} f_{j,i} =
%   \frac{1}{N_0} \sum_{i=1}^{N_0} f_{0,i}.\]

%%
% Define the subdivision operator that maps \(f_j\) to \(f_{j+1}\). 

subdivide = @(f,h)cconv( upsampling(f), h);

%%
% We use here the kernel
% \[ h = \frac{1}{8}[1, 4, 6, 4, 1]. \]
% It produced a cubic B-spline interpolation.

h = [1 4 6 4 1];
h = 2*h/sum(h(:));

%% 
% Initilize the subdivision with \(f_0\) at scale \(j=0\).

f = f0;

%%
% Perform one step.

f = subdivide(f,h);

%%
% Display the original and filtered discrete curves.

clf; hold on;
myplot(f, 'k.-');
myplot(f0, 'r.--');
myaxis(0);

%EXO
%% Perform several step of subdivision. Display the different curves.
Jmax = 3;
f = f0;
clf;
for j=0:Jmax
    f = subdivide(f,h);
    subplot(2,2,j+1);
    hold on;
    myplot(f, 'k.-');
    myplot(f0, 'r.--');
    myaxis(0);
end
%EXO

%%
% Under some restriction on the kernel \(h\), one can show that
% these discrete curves converges (e.g. in Hausdorff distance) toward a
% smooth limit curve \(f^\star : [0,1] \rightarrow \CC\). 

%%
% We do not details here sufficient condition to ensure convergence and
% smoothness of the limitting curve. The interested reader can have a look
% at <#biblio [DynLevin02]> for a review of theoritical guarantees for
% subdivision. 

%%
% The limit curve \(f^\star\) is a weighted average of the initial 
% points \(f_0 = (f_{0,i})_{i=0}^{N_0-1} \in \CC^{N_0}\) using a continuous 
% scaling function \(\phi : [0,1] \rightarrow \RR\)
% \[ f^\star(t) = \sum_{i=0}^{N_0-1} f_{0,i} \phi(t-i/N_0).  \]
% The continuous kernel \(\phi\) is a low-pass function which as a compact
% support of width \(K/N_0\). The control point \(f_{0,i}\) thus only
% influences the final curve \(f^\star\) around \(t=i/N_0\).

%%
% The scaling function \(\phi\) can be computed as the limit of the sub-division
% process \(f_j\) when starting from \(f_0 = \delta = [1,0,\ldots,0]\),
% which is the Dirac vector. 

%EXO
%% Compute the scaling function \(\phi\)
%% associated to the subdivision.
n = 6;
f = zeros(n,1); f(n/2+1) = 1;
for i=1:5
    f = subdivide(f,h);  
end
clf;
plot(linspace(-1/2,1/2,length(f)), f); axis([-1/2 1/2 -.01 max(f)*1.03]);
%EXO

%EXO
%% Test with different configurations of control points.
f0 = [0+0i; 1+0i; 1+1i; 0+1i];
f = f0;
clf;
for j=0:Jmax
    f = subdivide(f,h);
    subplot(2,2,j+1);
    hold on;
    myplot(f, 'k.-');
    myplot(f0, 'r.--');
    myaxis(.03);
end
%EXO

%% Quadratic B-splines
% We consider here the Chaikin "corner cutting"
% scheme <#biblio [Chaikin74]>.

%%
% For a weight \(w>1\), it corresponds to the following kernel:
% \[ h = \frac{1}{1+w}[1, w, w, 1]. \]
% The weight is a tension parameter that controls the properties of the
% interpolation.

hcc = @(w)[1 w w 1]/(1+w);

%%
% For \(w=3\), the scaling function \(\phi\) is a quadratic B-spline.


%EXO
%% Test the corner-cutting for \(w=3\).
h=hcc(3);
f = f0;
clf;
for j=0:Jmax
    f = subdivide(f,h);
    subplot(2,2,j+1);
    hold on;
    myplot(f, 'k.-');
    myplot(f0, 'r.--');
    myaxis(.03);
end
%EXO

%EXO
%% Test the corner-cutting for vaious values of \(w\).
wlist = [1 2 3 6];
subd = @(f,w)f;
for i=1:5
    subd = @(f,w)subdivide(subd(f,w),hcc(w));
end
lgd = {};
F = [];
for i=1:length(wlist);
    w = wlist(i);
    F(:,end+1) = subd(f0,w);
    lgd{i} = ['w=' num2str(w*16)];
end
clf;
hold on;
plot(F([1:end 1],:), 'LineWidth', 2);
myplot(f0, 'r.--');
legend(lgd);
myaxis(.03);
%EXO

%% Interpolating Subdivision
% Interpolating schemes keeps unchange the set of point at the previous
% level, and only smooth the position of the added points.

%%
% A subdivision is interpolating if the kernel satisfies
% \[ h(0)=1 \qandq \forall i \neq 0, \quad h(2i)=0. \]

%%
% We consider the four-point interpolation kernel
% proposed in <#biblio [DynLevGre87]>:
% \[ h = [-w, 0, 1/2+w, 1, 1/2+w, -w] \]
% where \(w>0\) is a tension parameter. 

h4pt = @(w)[-w, 0, 1/2+w, 1, 1/2+w, 0, -w];

%%
% One usually choose \(w=1/16\) wich corresponds to 
% cubic B-spline interpolation.
% It can be shown to produce \(C^1\) curves 
% for \( w \in [0, (\sqrt{5}-1)/8 \approx 0.154] \), see <#biblio [DynGreLev91]>.

%EXO
%% Perform the interpolating subdivision
%% for \(w=1/16\).
w = 1/16;
h = h4pt(w);
f = f0;
clf;
for j=0:Jmax
    f = subdivide(f,h);
    subplot(2,2,j+1);
    hold on;
    myplot(f, 'k.-');
    myplot(f0, 'r.--');    
    myaxis(.13);
    hold off;
end
%EXO

%EXO
%% Test the influence of \(w\).
wlist = [.5 1 1.5 2]/16;
subd = @(f,w)f;
for i=1:5
    subd = @(f,w)subdivide(subd(f,w),h4pt(w));
end
lgd = {};
F = [];
for i=1:length(wlist);
    w = wlist(i);
    F(:,end+1) = subd(f0,w);
    lgd{i} = ['w=' num2str(w*16) '/16'];
end
clf;
hold on;
plot(F([1:end 1],:), 'LineWidth', 2);
myplot(f0, 'r.--');
legend(lgd);
axis tight; axis off; axis equal;
%EXO

%EXO
%% Compare the result of the quadratic B-spline, cubic B-spline, 
%% and 4-points interpolating.
hh = [];
hh{end+1}=[1 3 3 1]/4;
hh{end+1} = [1 4 6 4 1]/8;
hh{end+1} = [-1, 0, 9, 1, 9, 0, -1]/16; 
hh{end}((end+1)/2)=1;
lgd = {'Quadratic', 'Cubic', 'Interpolating'};
col = {'r' 'g' 'b'};    
clf; hold on;
for k=1:length(hh)
    h = hh{k};
    f = f0;
    for j=0:7
        f = subdivide(f, h);
    end
    myplot(f, col{k});
end
myplot(f0, 'k.--');
axis tight; axis off;
legend(lgd);
%EXO


%%
% The 4-point scheme for \(w=1/16\) is generalized to \(2k\)-point schemes of
% Deslauriers-Dubuc
% <#biblio [DeslDub89]>. This corresponds to computing a polynomial
% interpolation of degree \(2k-1\), which generalizes the cubic interpolation.
% Using larger \(k\) leads to smoother interpolation, at the price of a
% larger interpolation kernel.

%%
% We give here the odd coefficients of the filters.

H = {   [0.5000 0.5000], ...
        [-0.0625, 0.5625, 0.5625, -0.0625], ...
        [0.0117, -0.0977, 0.5859, 0.5859, -0.0977, 0.0117], ...
        [-0.0024, 0.0239, -0.1196, 0.5981, 0.5981, -0.1196, 0.0239, -0.0024] };    
hdd = @(k)assign(assign(zeros(4*k-1,1),1:2:4*k-1,H{k}), 2*k, 1);

%EXO
%% Display the scaling function associated to these Deslauriers-Dubuc filters.
n = 8;
clf;
for k=1:4
    h = hdd(k);
    f = zeros(n,1); f(n/2+1) = 1;
    for i=1:5
        f = subdivide(f,h);
    end
    subplot(4,1,k);
    plot(linspace(-n/2,n/2,length(f)), f); 
    axis([-n/2 n/2 -.15 1.03]);
end
%EXO

    
%% Curve Approximation
% Given an input, complicated curve \(g : [0,1] \rightarrow \CC\), it is possible to approximate is by
% sampling the curve, and then subdividing it. It corresponds to a low pass
% filtering approximation.

%%
% Load an initial random curve, which is a high resolution curve \(g\).

options.bound = 'per';
n = 1024*2; 
sigma = n/8;
F = perform_blurring(randn(n,1),sigma,options) + 1i*perform_blurring(randn(n,1),sigma,options);
F = rescale(real(F),.01,.99) + 1i * rescale(imag(F),.01,.99);

%%
% Display it.

clf; myplot(F, 'k');
myaxis(0);

%%
% Load an interpolating subvision mask.

h = [-1, 0, 9, 1, 9, 0, -1]/16; 
h((end+1)/2)=1;

%EXO
%% Perform an approximation \(f\) of the curve using a uniform sampling with \(N_0=20\)
%% points.
p = 28;
t0 = (0:1/n:1-1/n)';
t = (0:1/p:1-1/p)';
f0 = interp1(t0,F,t);
f = f0;
Jmax = ceil(log2(n/p));
for j=0:Jmax
    f = subdivide(f, h);
end
clf; hold on;
myplot(F, 'k'); 
myplot(f0, 'k.');
myplot(f, 'r'); 
myaxis(0);
%EXO

%%
% To quantify the quality of the approximation, we use an averaged
% Hausdorff distance.
% The distance between two sets of points \(X\) and
% \(Y\) is
% \[ d(X,Y) = d_0(X,Y)+d_0(Y,X) \qwhereq
%       d_0(X,Y)^2 = \frac{1}{\abs{X}} \sum_{x \in X} \umin{y \in Y} \norm{x-y}^2. \]

%%
% Compute the pairwise distances matrix \(D_{i,j} = \norm{f_i-g_j}^2\) between points.

dist = @(f,g)abs( repmat(f, [1 length(g)]) - repmat(transpose(g), [length(f) 1]) );

%%
% Compute the Hausdorff distance.

hausdorff = @(f,g)sqrt( mean(min(dist(f,g)).^2) );
hausdorff = @(f,g)hausdorff(f,g) + hausdorff(g,f);

%EXO
%% Display the decay of the Hausdorff approximation error as the number \(N_0\) of
%% sampling points increases.
q = 15;
plist = round(linspace(8,80,q));
d = [];
for k=1:length(plist)
    p = plist(k);
    % sample
    t0 = (0:1/n:1-1/n)';
    t = (0:1/p:1-1/p)';
    f0 = interp1(t0,F,t);
    % interpolate
    f = f0;
    Jmax = ceil(log2(n/p));
    for j=0:Jmax
        f = subdivide(f, h);
    end
    % record distance    
    d(end+1) = hausdorff(F,f);
end
clf;
plot(plist, d, 'LineWidth', 2);
xlabel('N_0'); ylabel('d');
axis tight;
%EXO


%% 3-D Curve Subdivision
% It is possible to construct 3-D curves by subdivision.

%EXO
%% Perform curve subdivision in 3D space.
% control mesh.
f0 = rand(12,3);
Jmax = 4; ms = 20; lw = 1.5;
f = f0;
for j=0:Jmax
    f = cat(2, upsampling(f(:,1)), upsampling(f(:,2)), upsampling(f(:,3))  );
    f = cat(2, cconv(f(:,1),h), cconv(f(:,2),h), cconv(f(:,3),h) );    
end
clf;
subplot(1,2,1);
hold on;
hh = plot3([f(:,1);f(1,1)], [f(:,2);f(1,2)], [f(:,3);f(1,3)], 'k-');
set(hh, 'MarkerSize', ms);
set(hh, 'LineWidth', lw); 
hh = plot3([f0(:,1);f0(1,1)], [f0(:,2);f0(1,2)], [f0(:,3);f0(1,3)], 'r.--');
set(hh, 'LineWidth', lw);
axis('tight'); box('on'); view(3);  % axis('off');
subplot(1,2,2);
hold on;
hh = plot3([f(:,1);f(1,1)], [f(:,2);f(1,2)], [f(:,3);f(1,3)], 'k-');
set(hh, 'MarkerSize', ms);
set(hh, 'LineWidth', lw); 
hh = plot3([f0(:,1);f0(1,1)], [f0(:,2);f0(1,2)], [f0(:,3);f0(1,3)], 'r.--');
set(hh, 'LineWidth', lw);
axis('tight'); box('on');
view(70,25);
%EXO

%% Bibliography
% <html><a name="biblio"></a></html>


%%
% * [DynLevGre87] N. Dyn, D. Levin and J.A. Gregory, <http://dx.doi.org/10.1016/0167-8396(87)90001-X _A 4-point interpolatory subdivision scheme for curve design_>, Computer Aided Geometric Design, 4(4), Pages 257-268, 1987.
% * [Chaikin74] G. Chaikin, <http://dx.doi.org/10.1016/0146-664X(74)90028-8 _An algorithm for high speed curve generation_>. Computer Graphics and Image Processing, 3, 346-349, 1974.
% * [Riesen75] R. Riesenfeld, <http://dx.doi.org/10.1016/0146-664X(75)90017-9 _On Chaikin's algorithm_>. Computer Graphics and Image Processing 4, 3, 304-310, 1975.
% * [DeslDub89] G. Deslauriers and S. Dubuc. <http://dx.doi.org/10.1007/BF01889598 _Symmetric iterative interpolation processes_>. Constructive Approximation, 5(1):49-68, Dec. 1989.
% * [DynLevin02] N. Dyn and D. Levin. <http://dx.doi.org/10.1017/S0962492902000028 _Subdivision schemes in geometric modelling_>. Acta Numerica, 11:73-144, Jan. 2002.
% * [DynGreLev91] N. Dyn, J.A. Gregory, G. Levin, _Analysis of uniform binary subdivision schemes for curve design_, Constructive Approximation, 7(1), p. 127-147, 1991.
