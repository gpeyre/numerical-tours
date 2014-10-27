%% Curvelet Denoising
% This numerical tour explores the use of curvelets to perform image
% denoising.

perform_toolbox_installation('signal', 'general');

%CMT
rep = 'results/denoising_curvelets/';
if not(exist(rep))
    mkdir(rep);
end
%CMT

%% Curvelet Transform
% The curvelet tight frame was introduced by Candes and Donoho to enhance the wavelet 
% representation of geometric cartoon images. See for instance

%%
% E.J. Candès and D. L. Donoho,
% _New tight frames of curvelets and optimal representations of objects with piecewise-C2 singularities_,
% Comm. Pure Appl. Math., 57 219-266

%%
% In this tour, we use the discrete curvelet transform, detailed in 

%%
% E. J. Candès, L. Demanet, D. L. Donoho and L. Ying,
% _Fast discrete curvelet transforms_,
% Multiscale Model. Simul., 5  861-899

%%
% and it uses the <http://www.curvelet.org/ Curvelab> implementation.

%%
% Load an image.

n = 256;
name = 'lena';
M = rescale(load_image(name, n));


%%
% Parameters for the curvelet transform.

options.null = 0;
options.finest = 1;
options.nbscales = 4;
options.nbangles_coarse = 16;
options.is_real = 1;
options.n = n;

%%
% Perform the transform.

MW = perform_curvelet_transform(M, options);

%%
% Display the transform.

clf;
plot_curvelet(MW, options);

%%
% One can perform a non-linear approximation of the image by thresholding
% the curvelet coefficients.

T = .2;
MWT = perform_thresholding(MW, T, 'hard');

%%
% One recovers the approximated image by using the inverse curvelet
% transform.

M1 = perform_curvelet_transform(MWT, options);

%%
% Display the approximation. Note that is not the best non-linear
% approximation since the curvelet tight frame is not an orthogonal basis.

clf;
imageplot(M, 'Original', 1,2,1);
imageplot(clamp(M1), 'Approximated', 1,2,2);

%% Denoising with the Curvelet transform
% Curvelet thresholding is useful to perform denoising.

%%
% We load a sub-set of the lena image.

name = 'lena';
n = 128;
M0 = rescale(crop(load_image(name),n, [108 200]));
options.n = n;

%%
% Add some noise.

sigma = .05;
M = M0 + sigma*randn(n);

%EXO
%% Compute the best threshold to minimize the denoising error in curvelets.
%% Call |Mcurv| the optimal denoising.
Tlist = linspace(.8,1.2,15)*sigma;
MW = perform_curvelet_transform(M, options);
err = [];
for i=1:length(Tlist)
    MWT = perform_thresholding(MW, Tlist(i), 'hard');
    M1 = perform_curvelet_transform(MWT, options);
    err(end+1) = snr(M0,M1);
end
clf;
h = plot(Tlist/sigma,err);
set(h, 'LineWidth', 2);
axis('tight');
[tmp,i] = max(err);
MWT = perform_thresholding(MW, Tlist(i), 'hard');
Mcurv = perform_curvelet_transform(MWT, options);
%EXO

%%
% Display.

clf;
imageplot(clamp(M), 'Noisy', 1,2,1);
imageplot(clamp(Mcurv), ['Denoised, SNR=' num2str(snr(M0,Mcurv),3) 'dB'], 1,2,2);

%EXO
%% Perform cycle spinning to enhance the recovery error.
T = Tlist(i);
s = 4;
[dY,dX] = meshgrid(0:s-1,0:s-1);
Mcurv = zeros(n);
for i=1:s^2
    Ms = circshift(M, [dX(i) dY(i)]);
    MW = perform_curvelet_transform(Ms, options);
    MWT = perform_thresholding(MW, T, 'hard');
    Ms = perform_curvelet_transform(MWT, options);
    Ms = circshift(Ms, -[dX(i) dY(i)]);
    Mcurv = (1-1/i)*Mcurv + 1/i*Ms;
end
%EXO

%%
% Display.

clf;
imageplot(clamp(M), 'Noisy', 1,2,1);
imageplot(clamp(Mcurv), ['Denoised, SNR=' num2str(snr(M0,Mcurv),3)], 1,2,2);

%EXO
%% Compare with translation invariant hard thresholding.
options.ti = 1;
Jmin = 3;
T = 2.8*sigma;
MW = perform_wavelet_transf(M, Jmin, +1, options);
MWT = perform_thresholding(MW, T, 'hard');
Mwav = perform_wavelet_transf(MWT, Jmin, -1, options);
clf;
imageplot(clamp(Mwav), ['Wavelets, SNR=' num2str(snr(M0,Mwav),3) 'dB'], 1,2,1);
imageplot(clamp(Mcurv), ['Curvelet, SNR=' num2str(snr(M0,Mcurv),3) 'dB'], 1,2,2);
%EXO

%% Inverse Problem Regularization
% L1 sparsity in curvelets can be used to regularize inverse problems.

%EXO
%% Applies curvelet iterative thresholding to solve an inverse problem such
%% as inpainting, deconvolution or compressed sending.
%EXO
