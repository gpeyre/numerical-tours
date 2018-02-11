niter = 3000
lambda_list = linspace(.03, 0, niter)
err = []

for i in 1 : niter
    fTI = Psi(a)
    d = y - Phi(fTI, Omega)
    #step
    a = SoftThresh(a + tau*PsiS(Phi(d, Omega)) , lambda_list[i]*tau)
end
    
figure(figsize = (6, 6))
imageplot(clamP(fSpars), @sprintf("Sparsity inpainting TI, SNR = %.1f dB", snr(f0, fTI)))