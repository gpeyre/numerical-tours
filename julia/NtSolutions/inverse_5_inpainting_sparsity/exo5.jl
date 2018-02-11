lambda_list = linspace(1, 0, niter)
fHard = y

for i in 1 : niter
    fHard = Xi(HardThresh(PsiS(ProjC(fHard, Omega)), lambda_list[i]))
end

    
figure(figsize = (6, 6))
imageplot(clamP(fSpars), @sprintf("Inpainting hard thresh., SNR = %.1f dB", snr(f0, fHard)))