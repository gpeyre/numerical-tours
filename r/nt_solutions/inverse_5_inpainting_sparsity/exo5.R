lambda_list = seq(1, 0, length=niter)
fHard = y

for (i in 1:niter)
{
    fHard = Xi(HardThresh(PsiS(ProjC(fHard, Omega)), lambda_list[i]))
}
    
imageplot(fSpars, paste("Inpainting hard thresh., SNR =", round(snr(f0, fHard), 1), "dB"))