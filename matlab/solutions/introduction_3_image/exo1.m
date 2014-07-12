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
