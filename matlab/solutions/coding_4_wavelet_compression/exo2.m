% inverse transform
MT = perform_wavelet_transf(MWT,Jmin, -1);
MQ = perform_wavelet_transf(MWQ,Jmin, -1);
% error
eT = snr(M,MT);
eQ = snr(M,MQ);
% display
clf;
imageplot(MT, strcat(['Thresholding, SNR=' num2str(eT,2)]), 1,2,1)
imageplot(MT-MQ, strcat(['Thresholding - Approximating, SNR=+' num2str(eT-eQ,2)]), 1,2,2);
