options.ti = 1;
Jmin = 3;
T = 2.8*sigma;
MW = perform_wavelet_transf(M, Jmin, +1, options);
MWT = perform_thresholding(MW, T, 'hard');
Mwav = perform_wavelet_transf(MWT, Jmin, -1, options);
clf;
imageplot(clamp(Mwav), ['Wavelets, SNR=' num2str(snr(M0,Mwav),3) 'dB'], 1,2,1);
imageplot(clamp(Mcurv), ['Curvelet, SNR=' num2str(snr(M0,Mcurv),3) 'dB'], 1,2,2);
