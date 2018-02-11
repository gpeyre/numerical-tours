Jmin = 1;
M = .05*N^2;
KW = perform_wavelet_transf(K, Jmin, +1, options);
KWT = perform_thresholding(KW, M, 'largest');
clf;
plot_wavelet(abs(KWT)>0, Jmin, options);
