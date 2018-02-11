


ntests <- 4
tau_list <- seq(0.03, 0.25, length=ntests)
e0 <- -Inf
for (k in 1:ntests){
  tau <- tau_list[k]
  f1 <- NLmeans(tau)
  e <- snr(f0,f1)
  if (e>e0){
    fNL <- f1
    e0 <- e}
  imageplot(clamp(f1),paste("tau =", round(tau,4), ", SNR =", round(e,4), "dB"), c(2,2,k))
}
