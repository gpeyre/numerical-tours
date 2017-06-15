err_fft = collect(max(e, 1e-10) for e in (vecnorm(f, 2)^2 - cumsum(cR.^2)))

log10(err_fft./vecnorm(f)^2)[1:5000]
plot(log10(err_fft./vecnorm(f)^2), linewidth = 2)
title(L"$\log_{10}(\epsilon^2[M]/ ||f||^2)$")
xlim(1, n^2/50)
ylim(-2.35, 0)

show()