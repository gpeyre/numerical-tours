fL2 = PhiS(y);
clf;
e = snr(f0,fL2);
imageplot(clamp(fL2), ['SNR=' num2str(e,3) 'dB']);
