q = as.integer(sqrt(M))
fw = perform_wavelet_transf(f, Jmin + 1, +1)
fw1 = matrix(0, n, n)
fw1[1:q, 1:q] = fw[1:q, 1:q]
fM = perform_wavelet_transf(fw1, Jmin + 1, -1)
# Normalize to 0-1
fM = (fM-min(fM))/(max(fM)-min(fM))

imageplot(fM, paste("Linear, Wavelets, SNR = ",  round(snr(f, fM), 1), "dB"))