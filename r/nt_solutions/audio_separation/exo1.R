for (i in (1:p))
{
  Y[[i]] = perform_stft(y[,i],w,q,n)
  plot_spectogram(Y[[i]], paste("Source", i))
}