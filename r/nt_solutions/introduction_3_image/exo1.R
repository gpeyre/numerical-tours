# compute kernel h
t = linspace(-pi,pi,n)
h0 = matrix((cos(t) + 1) / 2 )
h = h0 %*% t(h0)
# compute FFT
Mf = fft(M[,] * h)
Lf = fftshift(log(abs(Mf) + 1e-1 ))
# display
imageplot(M[,] * h, 'Image', c(1,2,1))
imageplot(Lf, 'Fourier transform', c(1,2,2))