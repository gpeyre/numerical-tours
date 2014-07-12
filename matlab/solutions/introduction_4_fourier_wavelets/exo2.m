F = fft2(f);
a = sort(abs(F(:))); 
if a(1)<a(n0*n0)
	a = a(n0^2:-1:1);
end
T = a(M+1);
F = F .* (abs(F)>T);
fM = real( ifft2(F) );
% display
clf;
imageplot(clamp(fM), ['Non-linear, Fourier, SNR=' num2str(snr(f,fM), 4) 'dB']);
