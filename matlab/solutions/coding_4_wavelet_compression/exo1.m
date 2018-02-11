% thresholding approximation
MWT = perform_thresholding(MW, T, 'hard');
MWQ = perform_thresholding(MW, T, 'quantize');
% display
clf;
subplot(1,2,1);
plot_wavelet(MWT,Jmin);
title('Thresholded coefficients');
subplot(1,2,2);
plot_wavelet(MWT-MWQ,Jmin);
title('Thresholded - Quantized');
