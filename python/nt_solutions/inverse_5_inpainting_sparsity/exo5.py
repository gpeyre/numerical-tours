lambda_list = np.linspace(1, 0, niter)
fHard = y

for i in range(niter):
    fHard = Xi(HardThresh(PsiS(ProjC(fHard, Omega)), lambda_list[i]))

    
plt.figure(figsize=(6,6))
imageplot(clamp(fSpars), "Inpainting hard thresh., SNR = %.1f dB" %snr(f0, fHard))