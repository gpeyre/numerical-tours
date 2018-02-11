niter = 3000
lambda_list = seq(.03, 0, length=niter)
err = c()

for (i in 1:niter)
{
    fTI = Psi(a)
    d = y - Phi(fTI, Omega)
    #step
    a = SoftThresh(a + tau * PsiS(Phi(d, Omega)) , lambda_list[i] * tau) 
}

imageplot(clamp(fSpars), paste("Sparsity inpainting TI, SNR =", round(snr(f0, fTI), 1), "dB"))