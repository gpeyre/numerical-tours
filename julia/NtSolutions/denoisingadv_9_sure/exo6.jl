figure(figsize = (5, 5))
tmp = minimum(E)
val, i = findmin(E)
lambd = lambda_list[i]
hf = h(f, lambd)
imageplot(clamP(hf), @printf("SNR = %.1f dB", snr(f0, hf)))
