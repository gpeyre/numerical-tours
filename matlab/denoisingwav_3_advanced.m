%% Advanced Wavelet Thresholdings
% This numerical tour present some advanced method for denoising that makes
% use of some exoting 1D thresholding functions, that in some cases give
% better results than soft or hard thresholding.

perform_toolbox_installation('signal', 'general');

%% Generating a Noisy Image
% Here we use an additive Gaussian noise.

%%
% First we load an image.

n = 256;
name = 'hibiscus';
M0 = load_image(name,n);
M0 = rescale( sum(M0,3) );

%%
% Noise level.

sigma = .08;

%%
% Then we add some Gaussian noise to it.

M = M0 + sigma*randn(size(M0));

%%
% Compute a 2D orthogonal wavelet transform.

Jmin = 3;
MW = perform_wavelet_transf(M,Jmin,+1);

%% Semi-soft Thresholding
% Hard and soft thresholding are two specific non-linear diagonal
% estimator, but one can optimize the non-linearity to capture the
% distribution of wavelet coefficient of a class of images. 

%%
% Semi-soft thresholding is a familly of non-linearities that interpolates between soft and hard thresholding.
% It uses both a main threshold |T| and a secondary threshold |T1=mu*T|.
% When |mu=1|, the semi-soft thresholding performs a hard thresholding,
% whereas when |mu=infty|, it performs a soft thresholding.

T = 1; % threshold value
v = linspace(-5,5,1024);
clf;
hold('on');
plot(v, perform_thresholding(v,T,'hard'), 'b--');
plot(v, perform_thresholding(v,T,'soft'), 'r--');
plot(v, perform_thresholding(v,[T 2*T],'semisoft'), 'g');
plot(v, perform_thresholding(v,[T 4*T],'semisoft'), 'g:');
legend('hard', 'soft', 'semisoft, \mu=2', 'semisoft, \mu=4');
hold('off');

%EXO
%% Compute the denoising SNR for different values of |mu| and different
%% value of |T|. Important: to get good results, you should not threshold
%% the low frequency residual. 
mulist = linspace(1,12, 17);
nmu = length(mulist);
Tlist = linspace(.5,3.5,10)*sigma;
err = [];
for i=1:length(mulist)
    for j=1:length(Tlist)
        T = Tlist(j);
        MWT = perform_thresholding(MW,[T mulist(i)*T],'semisoft');
        MWT(1:2^Jmin,1:2^Jmin) = MW(1:2^Jmin,1:2^Jmin);
        MT = perform_wavelet_transf(MWT,Jmin,-1);
        err(i,j) = snr(M0,MT);
    end
end
clf;
imageplot(err); 
% Tlist/sigma, mulist, 
axis('tight');
set_label('T/\sigma', '\mu');
title('SNR of semi-soft thresholding');
axis('on');
%EXO

%%
% One can display, for each |mu|, the optimal SNR (obtained by testing many
% different |T|). For |mu=1|, one has the hard thresholding, which gives
% the worse SNR. The optimal SNR is atained here for |mu| approximately
% equal to 6.

err_mu = compute_min(err, 2);
clf;
plot(mulist, err_mu, '.-');
axis('tight');
set_label('\mu', 'SNR');

%% Soft and Stein Thresholding
% Another way to achieve a tradeoff between hard and soft thresholding is
% to use a soft-squared thresholding non-linearity, also named a Stein estimator.

%%
% We compute the thresholding curves.

T = 1; % threshold value
v = linspace(-4,4,1024);
% hard thresholding
v_hard = v .* (abs(v)>T);
% soft thresholding
v_soft = v .* max( 1-T./abs(v), 0 );
% Stein thresholding
v_stein = v .* max( 1-(T^2)./(v.^2), 0 );

%%
% We display the classical soft/hard thresholders.

clf;
hold('on');
plot(v, v_hard, 'b');
plot(v, v_soft, 'r');
plot(v, v_stein, 'k--');
hold('off');
legend('Hard', 'Soft', 'Stein');

%EXO
%% Compare the performance of Soft and Stein thresholders, by determining
%% the best threshold value.
Tlist = linspace(.5,4,20)*sigma;
err_hard = [];
err_soft = [];
err_stein = [];
for j=1:length(Tlist)
    T = Tlist(j);
    MWT = perform_thresholding(MW,T,'hard');
    MWT(1:2^Jmin,1:2^Jmin) = MW(1:2^Jmin,1:2^Jmin);
    MT = perform_wavelet_transf(MWT,Jmin,-1);
    err_hard(j) = snr(M0,MT);
    MWT = perform_thresholding(MW,T,'soft');
    MWT(1:2^Jmin,1:2^Jmin) = MW(1:2^Jmin,1:2^Jmin);
    MT = perform_wavelet_transf(MWT,Jmin,-1);
    err_soft(j) = snr(M0,MT);
    MWT = perform_thresholding(MW,T,'stein');
    MWT(1:2^Jmin,1:2^Jmin) = MW(1:2^Jmin,1:2^Jmin);
    MT = perform_wavelet_transf(MWT,Jmin,-1);
    err_stein(j) = snr(M0,MT);
end
clf;
hold('on');
plot(Tlist/sigma, err_hard, 'r');
plot(Tlist/sigma, err_soft, 'b');
plot(Tlist/sigma, err_stein, 'g');
hold('off');
axis('tight');
set_label('T/\sigma', 'SNR');
legend('Hard', 'Soft', 'Stein');
axis('on');
%EXO
