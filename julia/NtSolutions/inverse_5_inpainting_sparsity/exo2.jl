niter = 1000
lambda_list = linspace(.03, 0, niter)
err = []

for i in 1 : niter
    fSpars = SoftThreshPsi(ProjC(fSpars, Omega), lambda_list[i])
end
    
figure(figsize = (6, 6))
imageplot(clamP(fSpars), @sprintf("Sparsity inpainting, SNR = %.1f dB", snr(f0, fSpars)))