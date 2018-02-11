fw = NtToolBox.perform_wavelet_transf(f, Jmin, +1)
a = sort((abs(fw))[:])[end : -1 : 1] #sort a 1D copy of fw in descending order
T = a[M]
fw1 = fw .* (abs(fw) .> T)
fM = NtToolBox.perform_wavelet_transf(fw1, Jmin, -1)
figure(figsize = (5, 5))
imageplot(clamP(fM), @sprintf("Non-linear, Wavelets, SNR = %.1f dB", snr(f, fM)))