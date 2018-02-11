niter = 1000
lambda_list = seq(.03, 0, length=niter)
err = c()

for (i in 1:niter)
{
    fSpars = SoftThreshPsi(ProjC(fSpars, Omega), lambda_list[i])
}   

imageplot(clamp(fSpars), paste("Sparsity inpainting, SNR = ", round(snr(f0, fSpars), 1), "dB"))