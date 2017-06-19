Jmin = 1;
# forward transform
fW = perform_wavortho_transf(f,Jmin,+1,h);
# linear approximation
eta = 4;
fWLin = zeros(n,n);
fWLin[1:n/eta,1:n/eta] = fW[1:n/eta,1:n/eta];
# backward transform
fLin = perform_wavortho_transf(fWLin,Jmin,-1,h);
elin = snr(f,fLin);
# display
clf;
imageplot(f, "Original", 1,2,1); 
u = @sprintf("Linear, SNR=%.3f", elin);
imageplot(clamp(fLin), u, 1,2,2);
