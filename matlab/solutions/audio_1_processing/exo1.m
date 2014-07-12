% center for the Fourier analysis
t0 = n/4;
% width of the bump
sigma = 128;
% window
t = (1:n)';
h = exp( -(t-t0).^2 ./ (2*sigma^2) );
% fft
xh = x.*h;
tau = 1e-3;
xh = xh(abs(xh)>tau);
xf = fft(xh);
% display signal
clf; 
subplot(2,1,1);
plot(t, x); axis('tight');
set_graphic_sizes([], 20);
title('Signal x');
subplot(2,1,2);
plot(t, h); axis('tight');
set_graphic_sizes([], 20);
title('Window h');
% display FFTs
p = length(xf);
clf; 
subplot(2,1,1);
plot(xh);
set_graphic_sizes([], 20); axis('tight');
title('Windowed signal (zoom)');
subplot(2,1,2);
plot( -p/2+1:p/2, fftshift(abs(xf)) );
set_graphic_sizes([], 20); axis('tight');
title('Windowed FFT (zoom)');
