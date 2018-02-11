niter = 1000;
lambda_list = linspace(.03,0,niter);
err = [];
for i=1:niter
    fSpars = SoftThreshPsi( ProjC(fSpars,Omega), lambda_list(i) );    
end
clf;
imageplot(clamp(fSpars), ['Sparsity inpainting, SNR=' num2str(snr(f0,fSpars),3) 'dB']);
