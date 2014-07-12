f1 = perform_thresholding(fWSep,m,'largest');
for i=1:n
    f1(:,i) = perform_wavortho_transf(f1(:,i),Jmin,-1);
end
for i=1:n
    f1(i,:) = perform_wavortho_transf(f1(i,:)',Jmin,-1)';
end
% display
clf;
imageplot(clamp(Mnlin),strcat(['Isotropic, SNR=' num2str(snr(f,Mnlin),3) 'dB']), 1,2,1 );
imageplot(clamp(f1),strcat(['Separable, SNR=' num2str(snr(f,f1),3) 'dB']), 1,2,2 );
