M = .05*n^2;
KWT = perform_thresholding(KW, M, 'largest');
clf;
plot_wavelet(KWT~=0, Jmin, options);
title('Non zero coefficients.');
