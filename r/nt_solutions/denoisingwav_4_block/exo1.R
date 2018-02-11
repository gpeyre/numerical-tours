


Tlist <- seq(0.5, 2.5, length=30)*sigma
snr_stein <- vapply(Tlist, function(t){snr(f0, ThreshWav(f, t))}, double(1))

plot(Tlist/sigma, snr_stein, 'l', col='blue', xlab='T/sigma', ylab='SNR')

Tmax <- max(snr_stein)

fThresh = ThreshWav(f, T)
