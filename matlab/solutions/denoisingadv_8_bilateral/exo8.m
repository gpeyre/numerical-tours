Jmin = 3;
options.ti = 1;
fW = perform_wavelet_transf(f, Jmin, +1, options);
T = 2.8*mu;
fW = fW .* (abs(fW)>T);
fWav = perform_wavelet_transf(fW, Jmin, -1, options);
clf;
imageplot(clamp(fWav), ['Wavelets, SNR=' num2str(snr(f0,fWav),3) 'dB']);
