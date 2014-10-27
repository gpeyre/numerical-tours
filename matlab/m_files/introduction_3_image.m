%% Introduction to Image Processing
% This numerical tour explores some basic image processing tasks.


perform_toolbox_installation('signal', 'general');


%% Image Loading and Displaying
% Several functions are implemented to load and display images.

%%
% First we load an image.

% path to the images
name = 'lena';
n = 256;
M = load_image(name, []);
M = rescale(crop(M,n));

%%
% We can display it. It is possible to zoom on it, extract pixels, etc.

clf;
imageplot(M, 'Original', 1,2,1);
imageplot(crop(M,50), 'Zoom', 1,2,2);

%% Image Modification

%%
% An image is a 2D array, that can be modified as a matrix.

clf;
imageplot(-M, '-M', 1,2,1);
imageplot(M(n:-1:1,:), 'Flipped', 1,2,2);

%%
% Blurring is achieved by computing a convolution with a kernel.

% compute the low pass kernel
k = 9;
h = ones(k,k);
h = h/sum(h(:));
% compute the convolution
Mh = perform_convolution(M,h);
% display
clf;
imageplot(M, 'Image', 1,2,1);
imageplot(Mh, 'Blurred', 1,2,2);


%%
% Several differential and convolution operators are implemented.

G = grad(M);
clf;
imageplot(G(:,:,1), 'd/dx', 1,2,1);
imageplot(G(:,:,2), 'd/dy', 1,2,2);

%% Fourier Transform
% The 2D Fourier transform can be used to perform low pass approximation
% and interpolation (by zero padding).

%%
% Compute and display the Fourier transform (display over a log scale).
% The function |fftshift| is useful to put the 0 low frequency in the
% middle. After |fftshift|, the zero frequency is located at position
% (n/2+1,n/2+1).

Mf = fft2(M);
Lf = fftshift(log( abs(Mf)+1e-1 ));
clf;
imageplot(M, 'Image', 1,2,1);
imageplot(Lf, 'Fourier transform', 1,2,2);

%EXO
%% To avoid boundary artifacts and estimate really the frequency content of
%% the image (and not of the artifacts!), one needs to multiply |M| by a
%% smooth windowing function |h| and compute |fft2(M.*h)|. Use a sine
%% windowing function. Can you interpret the resulting filter ?
% compute kernel h
t = linspace(-pi(),pi(),n);
h = (cos(t)+1)/2;
h = h'*h;
% compute FFT
Mf = fft2(M.*h);
Lf = fftshift(log( abs(Mf)+1e-1 ));
% display
clf;
imageplot(M.*h, 'Image', 1,2,1);
imageplot(Lf, 'Fourier transform', 1,2,2);
%EXO

%EXO
%% Perform low pass filtering by removing the high frequencies of the
%% spectrum. What do you oberve ?
k = round(.8*n); k = round(k/2)*2; % even number
Mf = fft2(M);
Mf(n/2-k/2+2:n/2+k/2, n/2-k/2+2:n/2+k/2) = 0;
Mh = real( ifft2(Mf) );
% display
clf;
imageplot( crop(M), 'Image', 1,2,1);
imageplot(clamp( crop(Mh)), 'Low pass filtered', 1,2,2);
%EXO

%%
% It is possible to do image interpolating by adding high frequencies

p = 64;
n = p*4;
M = load_image('boat', 2*p); M = crop(M,p);
Mf = fftshift(fft2(M));
MF = zeros(n,n);
sel = n/2-p/2+1:n/2+p/2;
sel = sel;
MF(sel, sel) = Mf;
MF = fftshift(MF);
Mpad = real(ifft2(MF));
clf;
imageplot( crop(M), 'Image', 1,2,1);
imageplot( crop(Mpad), 'Interpolated', 1,2,2);

%%
% A better way to do interpolation is to use cubic-splines.
% It avoid ringing artifact because the spline kernel has a smaller support
% with less oscillations.

Mspline = image_resize(M,n,n);
clf;
imageplot( crop(Mpad), 'Fourier (sinc)', 1,2,1);
imageplot( crop(Mspline), 'Spline', 1,2,2);