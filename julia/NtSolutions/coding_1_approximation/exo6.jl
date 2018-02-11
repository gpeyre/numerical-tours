figure(figsize = (8, 8))
Mlist = [Int(round(n^2/100)), Int(round(n^2/20))]

for i in 1 : length(Mlist)
    M = Mlist[i]
    #perform hard thresholding 
    a = sort(abs(fC)[:])[end : -1 : 1] #sort a 1D copy of F in descending order
    T = a[M]
    fCT = fC.*(abs(fC) .> T)
    fM = plan_idct(fCT)*fCT
    imageplot(clamP(fM), @sprintf("M/ N = %.2f, SNR = %.1f dB", M/n^2,snr(f,fM)), [1, 2, i])
end