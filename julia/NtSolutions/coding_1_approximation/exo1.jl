figure(figsize = (8, 8))
Mlist = [round(n^2/100), round(n^2/20)]

for i in 1 : length(Mlist)
    M = Int(Mlist[i])
    #perform hard thresholding
    F = plan_fft(f)*f
    a = sort(reshape(abs(F), size(F)[1]*size(F)[2]))[end : -1 : 1] #sort a 1D copy of F in descending order
    T = a[M]
    fFT = F.*(abs(F) .> T)
    fM = real(plan_ifft(fFT)*fFT)
    imageplot(clamP(fM), @sprintf("M/N = %.2f, SNR = %.1f dB", M/n^2, snr(f, fM)), [1, 2, i])
end
