Tlist = linspace(.8,4.5,25)*sigma
err_soft = zeros([len(Tlist),1]) 
err_hard = zeros([len(Tlist),1])
for i in arange(0,len(Tlist)):
    aT = thresh_hard(a,Tlist[i])
    fWav = perform_wavortho_transf(aT,Jmin,-1,h)
    err_hard[i] = snr(f0,fWav)
    aT = thresh_soft(a,Tlist[i])
    aT[:2^Jmin:,:2^Jmin:] = a[:2^Jmin:,:2^Jmin:]
    fWav = perform_wavortho_transf(aT,Jmin,-1,h)
    err_soft[i] = snr(f0,fWav)
h1, = plot(Tlist/sigma,err_hard)
h2, = plot(Tlist/sigma,err_soft) 
axis('tight')
legend([h1,h2], ['Hard', 'Soft']);