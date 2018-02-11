




w <- array(rnorm(n0*n0), c(n0,n0)) / n0
w <- w - mean(w) + 1/n0**2

w_transf <- fft(w)
w_transf_extended <- array(0, c(n0,n0,1,dim(f0)[4]))
for (k in 1:dim(f0)[4]){
  w_transf_extended[,,,k] <- w_transf
}

FFT <- apply(f1, c(3,4), fft)
CONV <- array(0, dim(f1))
for (i in 1:n0){ for (j in 1:n0){ CONV[i,j,,] <- FFT[(i-1)*n0 + j,,] } }

CONV <- CONV*w_transf_extended


IFFT <- apply(CONV, c(3,4), fft, inverse=TRUE) / (n0*n0)
f <- array(0, dim(f1))
for (i in 1:n0){ for (j in 1:n0){ f[i,j,,] <- Re(IFFT[(i-1)*n0 + j,,]) } }




u <- array(1, c(n0,n0,1,3))
u[(round(n0/2) - round(n/2) + 1) : (round(n0/2) + round(n/2))
  ,(round(n0/2) - round(n/2) + 1) : (round(n0/2) + round(n/2)),,] <- f0
u[1,1,,] <- c(0,0,0)
u[2,1,,] <- c(1,1,1)
imageplot(as.cimg(clamp(u)), 'Input', c(1,2,1))



u <- f
u[1,1,,] <- c(0,0,0)
u[2,1,,] <- c(1,1,1)
imageplot(as.cimg(clamp(u)), 'Synthesized', c(1,2,2))

