%% Data Dependent Noise Models
% The simplest noise model is Gaussian additive noise, where the variance
% of the pixel value is independent of the mean (which is the value we look
% for). 

%%
% Unfortunately, most real-life data corresponds to noise model that are
% much more complicated, and in particular the variance of the noise often
% depends on the parameter of interest (for instance the mean).

perform_toolbox_installation('signal', 'general');

%% Poisson Noise
% A Poisson model assume that each pixel \(x\)
% of an image \(f(x)\) is drawn from
% a Poisson distribution of parameter \(\lambda=f_0(x)\), where
% \(f_0\) is the clean intensity image to recover.

%%
% \[ \PP(f(x)=k)=\lambda^k e^{-\lambda}/k! \]

%% 
% Display the Poisson distribution for several value of \(\lambda\).

lambda = [4 10 20];
[k,Lambda] = meshgrid(1:50, lambda);
P = Lambda.^k .* exp(-Lambda)./factorial(k);
h = plot(P'); axis('tight');
if using_matlab()
    set(h, 'LineWidth', 2);
end
legend('\lambda=2', '\lambda=10', '\lambda=20');
set_label('k', 'P(k)');

%%
% This model corresponds to a photon count, where \(\lambda\)
% is proportional to the number of photons that hits the receptor during
% the exposition time. This is useful to model medical imaging (confocal
% microscopy), TEP and SPECT tomography, and digital camera noises.

%%
% The goal of denoising is to retrieve the value of 
% \(f_0(x)=\lambda(x)\), the true
% image value, which corresponds to a clean image.

%%
% Note that \(\lambda\) is the mean of the distribution, but is also its variance, 
% so that the noise intensity perturbating the image pixel \(f(x)\)
% is proportional to \(f_0(x)\).

%%
% We load a clean, unquantized image.

n = 256;
name = 'lena';
f0u = rescale( load_image(name,n) );

%%
% Quantize the values to given \(\lambda\) range.

lmin = 1;
lmax = 40;
f0 = floor( rescale(f0u,lmin,lmax) );

%% 
% Generate the noisy image.

f = poissrnd(f0);

%%
% Display.

clf;
imageplot(f0, 'Intensity map f0', 1,2,1);
imageplot(f,  'Observed noisy image f', 1,2,2);

%%
% Display the difference, which shows that the noise level depends on the
% intensity. The noise is larger in bright areas.

clf;
imageplot(f0, 'Intensity map f0', 1,2,1);
imageplot(f-f0,  'f-f0', 1,2,2);

%EXO
%% Display noisy image contaminated by Poisson noise of varying range.
lmin = 1;
lmax = [5 10 50 100];
clf;
for i=1:length(lmax)
    f1 = poissrnd( floor( rescale(f0u,lmin,lmax(i)) ) );
    imageplot(f1, strcat(['\lambda_{max}=' num2str(lmax(i))]), 2,2,i );
end
%EXO

%%
% Parameters for the wavelet transform.

Jmin = 4;
options.ti = 1;


%%
% A translation invariance wavelet transform denoising computes 
% an estimate \(\tilde f\) of the clean image
% \(f_0\) as
% \[ \tilde f = \Psi^+ \circ S_T \circ \Psi f \]
% where \(\Psi\) is the wavelet transform and \(S_T\) the hard
% thresholding operator using a well chosen threshold \(T\). 

%EXO
%% Perform translation invariant wavelet hard thresholding directly on the
%% Poisson noisy image \(f\). Check for an optimal threshold that maximize
%% the SNR. Record the optimal result in |fPoisson|.
ntest = 30;
sigma = std(f(:)-f0(:));
Tlist = linspace(2,3.5,ntest)*sigma;
fW = perform_wavelet_transf(f,Jmin,+1,options);
err = [];
for i=1:ntest
    fWT = perform_thresholding(fW, Tlist(i), 'hard');
    f1 = perform_wavelet_transf(fWT,Jmin,-1,options);    
    err(i) = snr(f0,f1);
end
% display
clf;
h = plot(Tlist/sigma,err);
axis('tight');
if using_matlab()
    set(h, 'LineWidth', 2);
end
set_label('T/sigma', 'SNR');
% best result
[tmp,i] = max(err);
fWT = perform_thresholding(fW, Tlist(i), 'hard');
fPoisson = perform_wavelet_transf(fWT,Jmin,-1,options);
%EXO

%%
% Display.

clf;
imageplot(f0, 'Original image', 1,2,1);
imageplot(clamp(fPoisson,min(f0(:)),max(f0(:))), strcat(['Denoised, SNR=' num2str(snr(f0,fPoisson),4)]), 1,2,2);


%% Variance Stabilization Poisson Denoising
% A variance stabilization transform (VST)
% is a mapping \(\phi\) applied to a noisy image \(f\)
% so that the distribution of each pixel \(\phi(f(x))\) is
% approximately Gaussian.

%%
% The Anscombe variance stabilization transform (VST) is given by
% the non-linear mapping.

%%
% \[\phi(a) = 2 \sqrt{a+3/8}\]

%%
% It transforms a Poisson distribution \(P(\lambda)\) to a 
% distribution with approximate
% variance 1, whatever \(\lambda\) is, for \(\lambda\) large enough
% (say \(\lambda>20\)).

%%
% Other VST have been proposed, for instance the Freeman & Tukey VST, given
% by 

%%
% \[\phi(a) = \sqrt{a+1}+\sqrt{a}\]

%EXO
%% Display the estimated variance of a Poisson distribution for various 
%% \(\lambda\) (e.g. 10000 realization for each \(\lambda\)) and display
%% the variance of a stabilized distribution (here the green
%% curve corresponds to 'Freeman & Tukey' and the blue curve to 'Anscombe'.
lmax = 10;
ntest = 100000;
[V,U] = meshgrid(1:lmax,ones(ntest,1));
W = poissrnd(V);
Mstab1 = 2*sqrt(W+3/8);
Mstab2 = sqrt(W+1)+sqrt(W);
vstab = [];
vstab(1,:) = std(Mstab1,1).^2;
vstab(2,:) = std(Mstab2,1).^2;
v = std(W,1).^2;
% plot
clf;
%subplot(2,1,1);
%h = plot(v); axis('tight');
%if using_matlab()
%    set(h, 'LineWidth', 2);
%end
%title('Poisson variance');
%set_label('\lambda', 'Variance');
%subplot(2,1,2);
h = plot(vstab'); axis('tight');
if using_matlab()
    set(h, 'LineWidth', 2);
    legend('Anscombe', 'Freeman & Tukey');
end
set_label('\lambda', 'Variance');
%EXO

%% 
% To perform denosing, one applies the VST, performs
% wavelet thresholding as if the data was corrupted by an additive Gaussian
% noise of variance \(\sigma=1\), and then applies the inverse VST. 

%%
% This corresponds to computing an estimate \(\tilde f\) of the clean image
% \(f_0\) as
% \[ \tilde f = \phi^{-1} \pa{ \Psi^+ S_T \circ \Psi \circ \phi(f) } \]
% where \(\Psi\) is the wavelet transform and \(S_T\) the hard thresholding. 

%EXO
%% Perform translation invariance wavelet hard thresholding on the
%% variance stabilized image. Use for instance the Anscombe VST. 
%% Check for an optimal threshold that maximize
%% the SNR. Record the optimal result in |fVST|.
ntest = 30;
sigma = 1;
Tlist = linspace(2,3.5,ntest)*sigma;
fW = perform_wavelet_transf( 2*sqrt(f+3/8) ,Jmin,+1,options);
err = [];
for i=1:ntest
    fWT = perform_thresholding(fW, Tlist(i), 'hard');
    f1 = perform_wavelet_transf(fWT,Jmin,-1,options);    
    % undo VST
    f1 = (f1/2).^2 - 3/8;
    % record error
    err(i) = snr(f0,f1);
end
% display
clf;
h = plot(Tlist/sigma,err);
axis('tight');
if using_matlab()
    set(h, 'LineWidth', 2);
end
set_label('T/sigma', 'SNR');
% best result
[tmp,i] = max(err);
fWT = perform_thresholding(fW, Tlist(i), 'hard');
f1 = perform_wavelet_transf(fWT,Jmin,-1,options);
fVST = (f1/2).^2 - 3/8;
%EXO

%%
% Display.

clf;
imageplot(clamp(fPoisson,min(f0(:)),max(f0(:))), ...
        strcat(['Un-stabilized, SNR=' num2str(snr(f0,fPoisson),4)]), 1,2,1);
imageplot(clamp(fVST,min(f0(:)),max(f0(:))), ...
        strcat(['Stabilized, SNR=' num2str(snr(f0,fVST),4)]), 1,2,2);

%%
% There is no close form solution for the Freeman VST.
% For each denoised stabilized pixel value \(b \in \RR\), 
% one should use a Newton algorithm to solve the equation
% \(\phi(a)=b\) and recovered the denoised value \(a \in \RR\).
% This reads
% \[ a_{k+1} = a_k - \frac{\phi(a_k)-y}{\phi'(a_k)} \]

%EXO
%% Perform VST denoising using the Freeman VST.
%EXO

%% Multiplicative Noise
% A multiplicative noise corresponds to a noise model where the clean image
% \(f_0\) is multiplied by the noise \(W\) to obtained the noisy image \(f(x)=W(x) f_0(x)\).

%%
% This model is useful for data acquired with an active acquisition device,
% for instance SAR imaging and ultra-sound imaging.

%%
% The distribution of \(f(x)\) thus has mean \(f_0(x)\) and variance
% proportional to \(f_0(x) \sigma\) where \(\sigma\) is the variance of \(W\).

%%
% Load a clean image.

n = 256;
name = 'boat';
f0 = rescale( load_image(name,n), 1e-2,1 );

%% 
% A classical model for the multiplier noise \(W\), that is used for instance in
% SAR imaging, is to assume that the noise has a Gamma law 
% of mean 1 and variance parameterized by the noise level \(L\).
% \[ P(W=x) \sim x^{K-1} e^{ -x/\theta } \]

%%
% where the mean of \(P\) is \(s=K \theta\) and the variance is
% \(\sigma^2=K \theta^2\).

%%
% This corresponds to a acquisition device which is
% averaging \(K\) measures with 1-sided exponential distribution. 

K = 4;
sigma = 1/sqrt(K);

%%
% Generate the random multiplier.

W = gamrnd(1/sigma^2, sigma^2,n,n);

%% 
% Generate the noisy signal.

f = f0.*W;

%%
% Display.

clf;
imageplot(f0, 'Intensity map f0', 1,2,1);
imageplot(f,  'Observed noisy image f', 1,2,2);

%%
% Display the difference, which shows that the noise level depends on the
% intensity. The noise is larger in bright areas.

clf;
imageplot(f0, 'Intensity map f0', 1,2,1);
imageplot(f-f0,  'f-f0', 1,2,2);

%EXO
%% Generate several noisy images for several noise levels.
slist = [.1 .2 .3 .6];
clf;
for i=1:length(slist)
    Wu = gamrnd(1/slist(i)^2, slist(i)^2,n,n);
    imageplot(f0.*Wu, strcat(['\sigma=' num2str(slist(i))]), 2,2,i );
end
%EXO

%EXO
%% Perform translation invariance wavelet hard thresholding directly on the
%% noisy image \(f=f_0 W\). Check for an optimal threshold that maximize
%% the SNR. Record the optimal result in |fMult|.
ntest = 30;
sigma = std(f(:)-f0(:));
Tlist = linspace(2.8,4.5,ntest)*sigma;
fW = perform_wavelet_transf(f,Jmin,+1,options);
err = [];
for i=1:ntest
    fWT = perform_thresholding(fW, Tlist(i), 'hard');
    f1 = perform_wavelet_transf(fWT,Jmin,-1,options);    
    err(i) = snr(f0,f1);
end
% display
clf;
h = plot(Tlist/sigma,err);
axis('tight');
if using_matlab()
    set(h, 'LineWidth', 2);
end
set_label('T/sigma', 'SNR');
% best result
[tmp,i] = max(err);
fWT = perform_thresholding(fW, Tlist(i), 'hard');
fMult = perform_wavelet_transf(fWT,Jmin,-1,options);
%EXO

%%
% Display.

clf;
imageplot(f0, 'Original image', 1,2,1);
imageplot(clamp(fMult,min(f0(:)),max(f0(:))), strcat(['Denoised, SNR=' num2str(snr(f0,fMult),4)]), 1,2,2);

%% Variance Stabilization for Multiplicative Noise
% An approximate variance stabilization transform consist in taking the
% log. The \(\log(f)-a\) image is then equal to \(\log(f_0)\) contaminated by \(log(W)-a\) which is not
% too far from a centered Gaussian distribution if the variance of \(W\) is
% not too large. 

%%
% The value of \(a\) should be chosen as the mean value of the
% random variable \(\log(f)\) so that \(\log(W)-a\) has zero mean.
% There exists close form for the value of \(a\) as a function of \(\sigma=1/\sqrt{K}\),
% which we use here, where \(\psi\) is the polygamma function.

%%
% _Important:_ Scilab user should replace |psi| by |dlgamma|.

a = psi(K) - log(K);

%%
% Display stabililized image.

clf;
imageplot(f, 'f', 1,2,1);
imageplot(clamp(log(f)-a,-2,2), 'log(f)', 1,2,2);

%%
% Distribution of the noise multiplier and of the log.

clf;
subplot(2,1,1);
hist(W(:),100);
axis('tight');
subplot(2,1,2);
hist(log(W(:))-a,100);
axis('tight');

%EXO
%% Perform translation invariance wavelet hard thresholding on the
%% variance stabilized image using the log. 
%% Check for an optimal threshold that maximize
%% the SNR. Record the optimal result in |fVST|.
ntest = 30;
sigma = std(log(W(:)));
Tlist = linspace(2,3.5,ntest)*sigma;
fW = perform_wavelet_transf( log(f)-a,Jmin,+1,options);
err = [];
for i=1:ntest
    fWT = perform_thresholding(fW, Tlist(i), 'hard');
    f1 = perform_wavelet_transf(fWT,Jmin,-1,options);    
    % undo VST
    f1 = exp(f1);
    % record error
    err(i) = snr(f0,f1);
end
% display
clf;
h = plot(Tlist/sigma,err);
axis('tight');
if using_matlab()
    set(h, 'LineWidth', 2);
end
set_label('T/sigma', 'SNR');
% best result
[tmp,i] = max(err);
fWT = perform_thresholding(fW, Tlist(i), 'hard');
f1 = perform_wavelet_transf(fWT,Jmin,-1,options);
fVST = exp(f1);
%EXO

%%
% Display.

clf;
imageplot(clamp(fMult,min(f0(:)),max(f0(:))), strcat(['Un-stabilized, SNR=' num2str(snr(f0,fMult),4)]), 1,2,1);
imageplot(clamp(fVST,min(f0(:)),max(f0(:))), strcat(['Stabilized, SNR=' num2str(snr(f0,fVST),4)]), 1,2,2);
