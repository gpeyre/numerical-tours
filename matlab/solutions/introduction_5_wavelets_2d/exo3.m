M_list = round( n0^2 ./ [16 32 64 128] );
clf;
for i=1:length(M_list)
    M = M_list(i);
    f1 = PsiS( Thresh(fW,v(M+1)) );
    imageplot(clamp(f1), strcat(['M/N=1/' num2str(n0^2/M) ', SNR=' num2str(snr(f,f1),3) 'dB']), 2,2,i);
end
