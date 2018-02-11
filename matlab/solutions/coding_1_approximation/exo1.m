Mlist = round( [.01 .05]*n^2 );
clf;
for i=1:length(Mlist)
    M = Mlist(i);
    fFT = perform_thresholding(fF,M,'largest');
    fM = real( ifft2(fFT)*n );
    imageplot( clamp(fM), ['M/N=' num2str(M/n^2,2) ', SNR=' num2str(snr(f,fM),3) 'dB'], 1,2,i);
end
