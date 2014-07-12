niter = 3000;
lambda_list = linspace(.03,0,niter);
for i=1:niter
    fTI = Psi(a);
    d = y-Phi(fTI,Omega);
    % step 
    a = SoftThresh( a + tau * PsiS( Phi(d,Omega) ) , lambda_list(i)*tau ) ;
end
clf;
imageplot(clamp(fTI), ['Sparsity inpainting TI, SNR=' num2str(snr(f0,fTI),3) 'dB']);
