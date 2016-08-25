niter = 1000
lambda_list = np.linspace(.03, 0, niter)
err = []

for i in range(niter):
    fSpars = SoftThreshPsi(ProjC(fSpars, Omega), lambda_list[i])
    
plt.figure(figsize=(6,6))
imageplot(clamp(fSpars), "Sparsity inpainting, SNR = %.1f dB" %snr(f0, fSpars))