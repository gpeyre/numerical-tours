m = sqrt(M);
F = fftshift(fft2(f));
F1 = zeros(n0,n0);
sel = (n0/2-m/2:n0/2+m/2)+1;
F1(sel,sel) = F(sel,sel); 
fM = real( ifft2(fftshift(F1)) );
% display
clf;
imageplot(clamp(fM), ['Linear, Fourier, SNR=' num2str(snr(f,fM), 4) 'dB']);
