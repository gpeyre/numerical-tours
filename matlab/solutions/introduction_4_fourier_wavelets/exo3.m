q = sqrt(M);
fw = perform_wavelet_transf(f,Jmin,+1);
fw1 = zeros(n0,n0);
fw1(1:q,1:q) = fw(1:q,1:q);
fM = perform_wavelet_transf(fw1,Jmin,-1);
% display
clf;
imageplot(clamp(fM), ['Linear, Wavelets, SNR=' num2str(snr(f,fM), 4) 'dB']);
