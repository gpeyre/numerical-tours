f1 = zeros(n0,n0);
for i=1:tau_max*tau_max  
    fTrans = circshift(y,[X(i) Y(i)]);
    fTrans = PsiS( Thresh( Psi(fTrans) ,T) );
    fTrans = circshift(fTrans,-[X(i) Y(i)]);
    f1 = (i-1)/i*f1 + fTrans/i;
end
clf;
imageplot(clamp(y), 'Noisy image', 1,2,1);
imageplot(clamp(f1), strcat(['Denoising, SNR=' num2str(snr(f,f1),3) 'dB']), 1,2,2);
