clf;
for u=1:length(Mlist)
    M = Mlist(u);
    fLT = perform_thresholding(fL,M,'largest');
    % inverse
    fM = fLT;
    for i=1:n/w
        for j=1:n/w
            seli = (i-1)*w+1:i*w;
            selj = (j-1)*w+1:j*w;
            fM(seli,selj) = idct2( fM(seli,selj) );
        end
    end
    % display
    imageplot( clamp(fM), ['M/n^2=' num2str(M/n^2,2) ', SNR=' num2str(snr(f,fM),3) 'dB'], 1,2,u);
end
