



G <- grad(W)
G <- G[,,1] + 1i*G[,,2]
EvalG <- function(gamma){ bilinear_interpolate(G, Im(gamma), Re(gamma)) }
EvalW <- function(gamma){ bilinear_interpolate(W, Im(gamma), Re(gamma)) }
#
gamma <- gamma0
displist <- round(seq(1, niter, length=10))
k <- 1
imageplot(t(f))

for (i in 1:niter){
  n <- normal(gamma)
  g <- EvalW(gamma)*normalC(gamma) - n*dotp(EvalG(gamma), n)
  gamma <- resample( gamma + dt*g )
  if (i==displist[k]){
    lw <- 1
    if (i==1 || i==niter){
      lw <- 4
    }
    cplot(gamma, type="l", pch=1, lw=lw, lty=1, col="red", add=TRUE)
    k <- k+1
  }
}
