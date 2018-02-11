cR = sort(abs(fC)[:])[end : -1 : 1]
err_dct = collect(max(e, 1e-10) for e in (vecnorm(f, 2)^2 - cumsum(cR.^2)))

plot(log10(err_fft./vecnorm(f)^2), linewidth = 2, color = "red", label = "Fourier")
plot(log10(err_wav./vecnorm(f)^2), linewidth = 2, color = "green", label = "DCT")

title(L"$\log_{10}(\epsilon^2[M]/ ||f||^2)$")
xlim(1, n^2/50)
ylim(-2.35,0)
legend()

show()
