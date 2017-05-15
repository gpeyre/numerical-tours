t = linspace(-pi, pi, n)
h0 = (cos(t) + 1)./2
h = transpose(h0).*h0
# compute FFT
Mf = plan_fft(M.*h)*(M.*h)
Lf = fftshift(log( abs(Mf) + 1e-1 ))
# display
clf
imageplot(M.*h, "Image", [1, 2, 1])
imageplot(Lf, "Fourier transform", [1, 2, 2])
