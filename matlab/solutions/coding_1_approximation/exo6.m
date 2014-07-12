clf;
for i=1:length(Mlist)
    M = Mlist(i);
    fCT = perform_thresholding(fC,M,'largest');
    fM = idct2(fCT);
    imageplot( clamp(fM), ['M/N=' num2str(M/n^2,2) ', SNR=' num2str(snr(f,fM),3) 'dB'], 1,2,i);
end
