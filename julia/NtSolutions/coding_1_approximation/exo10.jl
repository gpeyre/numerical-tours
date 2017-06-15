figure(figsize = (8, 8))

for u in 1 : length(Mlist)
    M = Mlist[u]
    a = sort(abs(fL)[:])[end : -1 : 1] #sort a 1D copy of F in descending order
    T = a[M]
    fLT = fL.*(abs(fL) .> T)
    #fLT = perform_thresholding(fL, M, 'largest')
    fM = fLT
    
    for i in 1 : Base.div(n, w)
        for j in 1 : Base.div(n, w)
            seli = (i - 1)*w + 1 : i*w;
            selj = (j - 1)*w + 1 : j*w;
            fM[seli, selj] = plan_idct(fM[seli, selj])*fM[seli, selj];
        end
    end
    
    imageplot(clamP(fM), @sprintf("M/N = %.2f, SNR = %.1f dB", M/n^2, snr(f,fM)), [1, 2, u])
end
