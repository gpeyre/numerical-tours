niter = 300
Lambda_list = seq(1e-6, 0.01, length=20)

tau = 1.9 / ( 1 + max(Lambda_list) * 8 / epsilon)
fBest = y 
fTV = y

err = rep(0, length(Lambda_list))

for (it in  1:length(Lambda_list))
{
    Lambda = Lambda_list[it]

    for (i in 1:niter)
    {
        # Compute the gradient of the smoothed TV functional.
        Gr = grad(fTV)
        d = sqrt(epsilon**2 + (Gr**2)[,,1] + (Gr**2)[,,2])
        G = -div(Gr / array(rep(d, 2), dim=c(dim(d),2)))
        # step
        e = Phi(fTV, h) - y
        fTV = fTV - tau * ( Phi(e, h) + Lambda * G)
    }  
    err[it] = snr(f0, fTV)
    if (err[it] > snr(f0,fBest)){fBest = fTV}
}        

plot(Lambda_list, err, type="l", col=4, ylab="SNR", xlab="lambda")
fTV = fBest

imageplot(fBest, paste('TV deconvolution, SNR = ', round(snr(f0, fBest), 2), 'dB'))