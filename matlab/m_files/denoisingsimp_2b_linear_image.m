%% Linear Image Denoising
% This numerical tour introduces basic image denoising methods.

perform_toolbox_installation('signal', 'general');


%% Noisy Image Formation
% In these numerical tour, we simulate noisy acquisition by adding some
% white noise (each pixel is corrupted by adding an independant Gaussian
% variable). 

%%
% This is useful to test in an oracle maner the performance of our methods.

%%
% Size \(N = n \times n\) of the image.

n = 256;
N = n^2;

%% 
% We load a clean image \(x_0 \in \RR^N\).

name = 'hibiscus';
x0 = rescale( sum(load_image(name,n),3) );

%%
% Display the clean image.

clf;
imageplot(x0);

%%
% Variance of the noise.

sigma = .08;

%%
% We add some noise to it to obtain the noisy signal \(y = x_0 + w\).
% Here \(w\) is a realization of a Gaussian white noise of variance
% \(\si^2\).

y = x0 + sigma*randn(size(x0));

%%
% Display the clean and the noisy signals. 

clf;
imageplot(clamp(y));


%% Linear Image Denoising
% We consider a noising estimator \(x \in \RR^N\) of \(x_0\) that only
% depends on the observation \(y\). Mathematically speaking, it is thus a
% random vector that depends on the noise \(w\).

%%
% A translation invariant linear denoising is necessarely a convolution
% with a kernel \(h\)
% \[ x = x_0 \star h \]
% where the periodic convolution between two 2-D arrays is defined as
% \[ (a \star b)_i = \sum_j a(j) b(i-j). \]

%%
% It can be computed over the Fourier domain as
% \[ \forall \om, \quad \hat x(\om) = \hat x_0(\om) \hat h(\om). \]

cconv = @(a,b)real(ifft2( fft2(a).*fft2(b) ));

%%
% We use here a Gaussian fitler \(h\) parameterized by 
% the bandwith \(\mu\).

normalize = @(h)h/sum(h(:));
t = [0:n/2-1, -n/2:-1]';
[Y,X] = meshgrid(t,t);
h = @(mu)normalize( exp( -(X.^2+Y.^2)/(2*mu^2) ) );

%%
% Display the filter \(h\) and its Fourier transform.

mu = 10;
clf;
subplot(2,1,1);
imageplot( fftshift( h(mu) )  ); axis('tight');
title('h');
subplot(2,1,2);
imageplot( fftshift( real(fft2(h(mu))) )  ); axis('tight');
title('fft2(h)');

%%
% Shortcut for the convolution with \(h\).

denoise = @(x,mu)cconv(h(mu), x);

%%
% Display a denoised signal.

clf;
imageplot( denoise(y,mu) );

%EXO
%% Display a denoised signal for several values of \(\mu\).
mu_list = linspace(.5,6,6);
clf;
for i=1:length(mu_list)
    mu = mu_list(i);
    subplot(2,3,i);
    imageplot( denoise(y,mu) );
    title(['\mu=' num2str(mu)]);
end
%EXO

%EXO
%% Display the evolution of the oracle denoising error 
%% \( \norm{y-x_0} \) as a function of \(\mu\).
%% Set \(\mu\) to the value of the optimal parameter.
mulist = linspace(.1,3.5,31);
err = arrayfun(@(mu)norm(x0-denoise(y,mu), 'fro' ), mulist);
clf;
h = plot(mulist,err); axis('tight');
set_label('\mu', '|y-x_0|');
% retrieve the best denoising result
[~,i] = min(err);
mu = mulist(i);
%EXO

%%
% Display the results.

clf;
imageplot( denoise(y,mu) );

%% Wiener Filtering
% We suppose here that \(x_0\) is a realization of a random vector \(x_0\), 
% whose distribution is Gaussian with a stationary covariance \(c\), 
% and we denote \(P_{X_0}(\om) = \hat c(\om)\) the power-spectrum of 
% \(x_0\). 

%%
% Recall that \(w\) is a realization of a random vector \(W\)
% distributed according to \(\Nn(0,\si^2 \text{Id})\).

%%
% The (oracle) optimal filter minimizes the risk 
% \[ R(h) = \EE_{W,X_0}( \norm{ X_0 - h \star (X_0 + W) }^2 ). \]

%%
% One can show that the solution of this problem, the so-called Wiener filter, 
% is defined as 
% \[ \forall \om, \quad \hat h(\om) = \frac{ P_{X_0}(\om) }{ P_{X_0}(\om) + \si^2 }. \]

%%
% We estimate \( P_{X_0} \) using the periodogram associated to the
% realization \(x_0\), i.e.
% \[ P_{X_0} \approx \frac{1}{N} \abs{\hat x_0}^2. \]

P = 1/N * abs(fft2(x0)).^2;

%%
% Compute the approximate Wiener filter.

h_w = real( ifft2( P ./ ( P + sigma^2 ) ) );

%%
% Note that this is a theoretical filter, because in practice one does not
% have access to \(x_0\).

%%
% Display it.

clf;
imageplot(crop(fftshift(h_w),n/8)); axis tight;

%%
% Display the denoising result.

clf;
imageplot(cconv(y,h_w));

%%
% Note that this denoising is not very efficient, because the hypothesis of
% stationarity of \(X_0\) is not realistic for such piecewise-regular
% signal. 