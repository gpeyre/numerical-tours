


tlist <- seq(.5, 2, length=30)*sigma

snr_stein <- c()
for (t in tlist){
  snr_stein <- c(snr_stein, snr(f0, ThreshWav(f,t)))
}

plot(tlist/sigma, snr_stein, 'l', col='blue', xlab='T/sigma', ylab='SNR')

Tmax <- max(snr_stein)
fBlock <- ThreshWav(f, T)

