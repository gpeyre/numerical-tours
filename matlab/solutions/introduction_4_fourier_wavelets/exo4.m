fw = perform_wavelet_transf(f,Jmin,+1);
a = sort(abs(fw(:))); 
if a(1)<a(n0*n0)
	a = a(n0^2:-1:1);
end
T = a(M+1);
fw1 = fw .* (abs(fw)>T);
fM = perform_wavelet_transf(fw1,Jmin,-1);
% display
clf;
imageplot(clamp(fM), ['Non-linear, wavelets, SNR=' num2str(snr(f,fM), 4) 'dB']);
