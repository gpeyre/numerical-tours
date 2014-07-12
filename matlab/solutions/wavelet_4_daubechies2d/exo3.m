fLin = perform_wavortho_transf(fWLin,0,-1);
elin = snr(f,fLin);
clf;
imageplot(f, 'Original', 1,2,1); 
imageplot(clamp(fLin), strcat(['Linear, SNR=' num2str(elin)]), 1,2,2);
