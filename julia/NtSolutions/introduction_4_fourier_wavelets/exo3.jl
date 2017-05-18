q = Int(sqrt(M))
fw = perform_wavelet_transf(f, Jmin, +1)
fw1 = zeros(n0, n0)
fw1[1:q, 1:q] = fw[1:q, 1:q]
fM = perform_wavelet_transf(fw1, Jmin, -1)
figure(figsize = (5, 5))
imageplot(clamP(fM),  @sprintf("Linear, Wavelets, SNR = %.1f dB", snr(f, fM)))
