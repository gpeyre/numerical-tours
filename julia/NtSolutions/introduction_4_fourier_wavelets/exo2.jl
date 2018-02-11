F = plan_fft(f)*f
a = sort(reshape(abs(F), size(F)[1]*size(F)[2]))[end : -1 : 1] #sort a 1D copy of F in descending order
T = a[M]
FT = F.*(abs(F) .> T)
fM = real(plan_ifft(FT)*FT)
figure(figsize = (5,5))
imageplot(clamP(fM), @sprintf("Non-linear, Fourier, SNR = %.1f dB", snr(f, fM)))
