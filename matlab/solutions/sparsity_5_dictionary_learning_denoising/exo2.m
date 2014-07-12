Jmin = 3;
options.ti = 1;
MW = perform_wavelet_transf(f, Jmin, +1, options);
MW = perform_thresholding(MW, 2.8*sigma, 'hard');
Mwav = perform_wavelet_transf(MW, Jmin, -1, options);
clf;
imageplot(clamp(Mwav), ['Wavelets, SNR=' num2str(snr(f0,Mwav),4) 'dB']);
