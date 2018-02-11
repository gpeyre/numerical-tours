niter = 500;
lambda_list = linspace(1,0,niter);
fHard = y; 
for i=1:niter
    fHard = Xi( HardThresh(PsiS(ProjC(fHard,Omega)), lambda_list(i)) );
end
clf;
imageplot(clamp(fHard), ['Inpainting hard thresh., SNR=' num2str(snr(f0,fHard),3) 'dB']);
