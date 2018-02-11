options.ti = 1;
Jmin = 4;
MW = perform_wavelet_transf(M, Jmin, +1, options);
MWT = perform_thresholding(MW, 3*sigma, 'hard');
Mwav = perform_wavelet_transf(MWT, Jmin, -1, options);
clf;
imageplot(clamp(Mnl1), ['NL-means, SNR=' num2str(snr(M0,Mnl1)) 'dB'], 1,2,1);
imageplot(clamp(Mwav), ['Wavelets, SNR=' num2str(snr(M0,Mwav)) 'dB'], 1,2,2);
