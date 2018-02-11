tlist = linspace(.5, 2, 20).*sigma

snr_stein = collect(snr(f0, ThreshWav(f,t)) for t in Tlist)

plot(collect(T/sigma for T in Tlist), snr_stein, linewidth = 2)
xlabel(L"$T/sigma$")
ylabel("SNR")

show()

Tmax = maximum(snr_stein)
fTI = ThreshWav(f, T)