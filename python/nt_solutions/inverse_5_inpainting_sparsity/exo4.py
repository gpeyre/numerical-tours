niter = 3000
lambda_list = np.linspace(.03, 0, niter)
err = []

for i in range(niter):
    fTI = Psi(a)
    d = y-Phi(fTI, Omega)
    #step
    a = SoftThresh(a + tau*PsiS(Phi(d, Omega)) , lambda_list[i]*tau) 
    
plt.figure(figsize=(6,6))
imageplot(clamp(fSpars), "Sparsity inpainting TI, SNR = %.1f dB" %snr(f0, fTI))