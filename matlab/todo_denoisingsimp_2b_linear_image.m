%% Linear Sigal Denoising
% This numerical tour introduces basic image denoising methods.

perform_toolbox_installation('signal', 'general');

%% Noisy Signal Formation
% In these numerical tour, we simulate noisy acquisition by adding some
% white noise (each pixel is corrupted by adding an independant Gaussian
% variable). 

%%
% This is useful to test in an oracle maner the performance of our methods.

%%
% Length \(N\) of the signal.

N = 1024;

%% 
% We load a clean signal \(x_0 \in \RR^N\).

name = 'piece-regular';
x0 = rescale( load_signal(name,N) );

%%
% We add some noise to it to obtain the noisy signal \(x = x_0 + w\).
% Here \(w\) is a realization of a Gaussian white noise of variance
% \(\si^2\).

sigma = .04; % noise level
x = x0 + sigma*randn(size(x0));
clf;
subplot(2,1,1);
plot(x0); axis([1 N -.05 1.05]);
subplot(2,1,2);
plot(x); axis([1 N -.05 1.05]);

%%
% We load an image.

N = 256;
name = 'hibiscus';
M0 = load_image(name,N);
M0 = rescale( sum(M0,3) );

%%
% Then we add some gaussian noise to it.

sigma = .08; % noise level
M = M0 + sigma*randn(size(M0));
clf;
imageplot(M0, 'Original', 1,2,1);
imageplot(clamp(M), 'Noisy', 1,2,2);

%% Linear Signal Denoising
% A translation invariant linear denoising is necessarely a convolution
% with a kernel |h|. It correspond to a linear diagonal operation over the
% Fourier domain that multiply each noisy Fourier coefficient by the
% Fourier transform of |h|.

%%
% In practice, one uses a Gaussian fitler |h|, and the only parameter
% is the width (variance) of the filter. 

% width of the filter
mu = 4;
% compute a Gaussian filter of width mu
t = (-length(x)/2:length(x)/2-1)';
h = exp( -(t.^2)/(2*mu^2) );
h = h/sum(h(:));

%%
% The Fourier transform of a Gaussian discrete filter is nearly a Gaussian
% whose width is proportional to |1/mu|.

% Fourier transform of the (centered) filter
hf = real( fft(fftshift(h)) );
hf = fftshift(hf);
% display
clf;
subplot(2,1,1);
plot( t,h  ); axis('tight');
title('Filter h');
subplot(2,1,2);
plot( t,hf  ); axis('tight');
title('Fourier transform');


%% 
% Since we use periodic boundary condition, the convolution of |x| with |h|
% can be computer over the Fourier domain.

% Fourier coefficients of the noisy signal
xf = fft(x);
% Fourier coefficients of the denoised signal
xhf = xf .* fft(fftshift(h));
% Denoised signal
xh = real( ifft(xhf) );

%%
% We display the denoised signal. 
% Although most of the noise is removed, the singularity have been blurred.

clf;
subplot(2,1,1);
plot( t,x  ); axis('tight');
title('Noisy');
subplot(2,1,2);
plot( t,xh  ); axis('tight');
title('Denoised');


%%
% We display the noisy and denoised Fourier coefficients.
% One can see that the denoising remove the high frequency coefficients.

% log of Fourier transforms
epsilon = 1e-10;
L0 = log10(epsilon+abs(fftshift(fft(x0))));
L = log10(epsilon+abs(fftshift(xf)));
Lh = log10(epsilon+abs(fftshift(xhf)));
% display Fourier transforms
clf;
subplot(2,1,1);
plot( t, L, '-' );
axis([-length(x)/2 length(x)/2 -4 max(L)]);
title('log of noisy Fourier coefs.');
subplot(2,1,2);
plot( t, Lh, '-' );
axis([-length(x)/2 length(x)/2 -4 max(L)]);
title('log of denoised Fourier coefs.');

%% 
% It is non-trivial to select the width 
% parameter |mu| to minimize the denoising error.
% It should account for both the variance of the noise and
% the power spectrum of the image. 

%%
% We display the blurring for an increasing value of |mu|.

mulist = linspace(.5,4,4);
clf;
for i=1:length(mulist)
    mu = mulist(i);
    % compute the filter   
    h = exp( -(t.^2)/(2*mu^2) );
    h = h/sum(h(:));
    % perform the blurring
    xh = real( ifft( fft(x) .* fft(fftshift(h)) ));
    subplot( 4,1,i );
    plot(t, clamp(xh) ); axis('tight');
    title(strcat(['\mu=' num2str(mu)]));
end

%EXO
%% Try for various Gaussian variance |mu| to compute the denoising |xh|.
%% Compute, in an oracle manner, the best variance |muopt| by computing the
%% residual error |snr(x0,xh)|.
% display blurring for various mu
% compute the error for many mu
mulist = linspace(.1,3.5,31);
err =  [];
for i=1:length(mulist)
    mu = mulist(i);
    % compute the filter   
    h = exp( -(t.^2)/(2*mu^2) );
    h = h/sum(h(:));
    % perform blurring
    xh = real( ifft( fft(x) .* fft(fftshift(h)) ));
    err(i) = snr(x0,xh);
end
clf;
plot(mulist,err, '.-'); axis('tight');
set_label('\mu', 'SNR');
% retrieve the best denoising result
[snr_opt,I] = max(err);
muopt = mulist(I);
disp( strcat(['The optimal smoothing width is ' num2str(muopt) ' pixels, SNR=' num2str(snr_opt) 'dB.']) );
%EXO


%%
% Display the results.

% compute the optimal filter
h = exp( -(t.^2)/(2*muopt^2) );
h = h/sum(h(:));
% perform blurring
xh = real( ifft( fft(x) .* fft(fftshift(h)) ));
% display
clf;
subplot(2,1,1);
plot(t, clamp(x)); axis('tight');
title('Noisy');
subplot(2,1,2);
plot(t, clamp(xh)); axis('tight');
title('Denoised');


%% Linear Image Denoising
% We denoise similarely a 2D image using a 2D Gaussian filter whose width
% |mu| is optimized to match the noise level and the regularity of the
% signal.


%%
% We use a simple gaussian blur to denoise an image.

% we use cyclic boundary condition since it is quite faster
options.bound = 'per';
% number of pixel of the filter
mu = 10;
Mh = perform_blurring(M,mu,options);
clf;
imageplot(clamp(M), 'Noisy', 1,2,1);
imageplot(clamp(Mh), 'Blurred', 1,2,2);

%%
% We display the blurring for an increasing value of |mu|.

mulist = linspace(3,15,6);
clf;
for i=1:length(mulist)
    mu = mulist(i);
    Mh = perform_blurring(M,mu,options);
    imageplot(clamp(Mh), strcat(['\mu=' num2str(mu)]), 2,3,i);    
end

%EXO
%% Try for various Gaussian variance to compute the denoising |Mh|.
%% Compute, in an oracle manner, the best variance |muopt| by computing the
%% residual error |snr(M0,Mh)|.
% now compute the error for many mu
mulist = linspace(.3,6,31);
err =  [];
for i=1:length(mulist)
    mu = mulist(i);
    Mh = perform_blurring(M,mu,options);
    err(i) = snr(M0,Mh);
end
clf;
plot(mulist,err, '.-'); axis('tight');
set_label('\mu', 'SNR');
% retrieve the best denoising result
[snr_opt,I] = max(err);
muopt = mulist(I);
disp( strcat(['The optimal smoothing width is ' num2str(muopt) ' pixels, SNR=' num2str(snr_opt) 'dB.']) );
%EXO

%%
% Display the results

% optimal filter
Mgauss = perform_blurring(M,muopt,options);
% display
clf;
imageplot(M, strcat(['Noisy, SNR=' num2str(snr(M0,M)) 'dB']), 1,2,1);
imageplot(Mgauss, strcat(['Gaussian denoise, SNR=' num2str(snr(M0,Mgauss)) 'dB']), 1,2,2);

%% Wiener filtering
% In a probabilistic setting, for translation invariant signal
% distributions, the Wiener filtering is the optimal filtering.

%%
% Perform the wiener filtering

% FFT-based wiener filtering (using the oracle fourier coefficients)
x0f = fft2(M0);
Pxf = abs(x0f).^2; % power spectra
Hf = Pxf./(Pxf + N*sigma^2); % filter fourier transform
% compute convolution
xf = fft2(M);
Mwien = real( ifft2(xf.*Hf) );
Hwien = real( fftshift( ifft2(Hf) ) );

%% 
% display the filter

k = 5;
clf;
imageplot(Hwien(N/2-k+2:N/2+k,N/2-k+2:N/2+k), 'Wiener filter (zoom)');

%%
% display the result

% display
clf;
imageplot( clamp(Mgauss), strcat(['Gaussian denoise, SNR=' num2str(snr(M0,Mgauss)) 'dB']), 1,2,1);
imageplot( clamp(Mwien), strcat(['Wiener denoise, SNR=' num2str(snr(M0,Mwien)) 'dB']), 1,2,2);