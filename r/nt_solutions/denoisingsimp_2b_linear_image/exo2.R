

mulist <- seq(0.1,3.5, length=31)
err <- matrix( 0, c(length(mulist),1) )
for ( i in (1:length(mulist)) ){
  mu <- mulist[i]
  err[i] <- norm( as.matrix(x0) - denoise(as.matrix(y), mu) )
}

plot(mulist, err, 'l')
# retrieve the best denoising result
i <- which.min(err)
mu <- mulist[i]