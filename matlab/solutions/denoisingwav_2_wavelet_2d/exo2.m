fTI = zeros(n,n);
T = 3*sigma;
for i=1:m^2
    fS = circshift(f,[dX(i) dY(i)]);
    a = perform_wavelet_transf(fS,Jmin,1,options);
    aT = perform_thresholding(a,T,'hard');    
    fS = perform_wavelet_transf(aT,Jmin,-1,options);
    fS = circshift(fS,-[dX(i) dY(i)]);
    fTI = (i-1)/i*fTI + 1/i*fS;
end
clf;
imageplot(clamp(fHard), strcat(['Hard denoising, SNR=' num2str(snr(f0,fHard),3)]), 1,2,1);
imageplot(clamp(fTI), strcat(['Cycle spinning denoising, SNR=' num2str(snr(f0,fTI),3)]), 1,2,2);
