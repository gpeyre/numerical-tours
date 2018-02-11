q = Int(sqrt(M))
F = fftshift(plan_fft(f)*f)
Sel = zeros(n0, n0)

Sel[div(n0, 2) - div(q, 2) : div(n0, 2) + div(q, 2), div(n0, 2) - div(q, 2) : div(n0, 2) + div(q, 2)] = 1
F_zeros = F.*Sel

fM = real(plan_fft(fftshift(F_zeros))fftshift(F_zeros))
figure(figsize = (5,5))
imageplot(clamP(fM), @sprintf("Linear, Fourier, SNR = %.1f dB", snr(f, fM)))
