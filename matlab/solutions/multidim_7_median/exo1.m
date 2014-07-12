Mindep = zeros(n,n,3);
for i=1:3
    Mindep(:,:,i) = perform_median_filtering(M(:,:,i),k);    
end
pindep = snr(M0,Mindep);
clf;
imageplot(clamp(M), strcat(['Noisy, SNR=' num2str(pnoisy)]), 1,2,1 );
imageplot(clamp(Mindep), strcat(['Denoised, SNR=' num2str(pindep)]), 1,2,2 );
