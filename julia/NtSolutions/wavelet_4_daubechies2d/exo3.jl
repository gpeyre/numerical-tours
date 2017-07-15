Jmin = 1;
# forward transform
fW = perform_wavortho_transf(f, Jmin, +1, h);
# linear approximation
eta = 4;
fWLin = zeros(n, n);
fWLin[1 : Int(n/eta), 1 : Int(n/eta)] = fW[1 : Int(n/eta), 1 : Int(n/eta)];
# backward transform
fLin = perform_wavortho_transf(fWLin, Jmin, -1, h);
elin = snr(f, fLin);
# display
clf;
imageplot(f, "Original", [1, 2, 1]);
u = @sprintf("Linear, SNR=%.3f", elin);
imageplot(clamP(fLin), u, [1, 2, 2]);
