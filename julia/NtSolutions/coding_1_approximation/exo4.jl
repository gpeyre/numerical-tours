figure(figsize = (8, 8))
Mlist = [Int(round(n^2/100)), Int(round(n^2/20))]

for i in 1 : length(Mlist)
    M = Mlist[i]
    a = sort(abs(fW)[:])[end:-1:1] #sort a 1D copy of F in descending order
    T = a[M]
    fWT = fW.*(abs(fW) .> T)
    fM = real(perform_wavortho_transf(fWT, Jmin, -1, h))
    imageplot(clamP(fM), @sprintf("M/N = %.2f, SNR = %.1f dB", M/n^2, snr(f,fM)),[1, 2, i])
end