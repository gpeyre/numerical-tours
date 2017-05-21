Tlist = linspace(.5, 2.5, 30).*sigma
snr_stein = collect(snr(f0, ThreshWav(f, T)) for T in Tlist)

plot(collect(T/sigma for T in Tlist), snr_stein, linewidth = 2)
xlabel(L"$T/\sigma$")
ylabel("SNR")

show()

Tmax = maximum(snr_stein)

fThresh = ThreshWav(f, T)
