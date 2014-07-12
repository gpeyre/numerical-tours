clf;
for i=1:length(Mlist)
    M = Mlist(i);
    fWT = perform_thresholding(fW,M,'largest');
    fM = perform_wavortho_transf(fWT,Jmin,-1, options);
    imageplot( clamp(fM), ['M/N=' num2str(M/n^2,2) ', SNR=' num2str(snr(f,fM),3) 'dB'], 1,2,i);
end
