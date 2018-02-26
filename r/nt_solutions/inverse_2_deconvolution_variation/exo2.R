Lambda_list = seq(2*1e-3, 0.7, length=40)
nbiter = 4000

err = rep(0, length(Lambda_list))

for (i in 1:length(Lambda_list))
{
    Lambda = Lambda_list[i]
    fSob = Re(fft(yF * hF / (abs(hF)**2 + Lambda*S), inverse=TRUE)) / length(yF * hF / (abs(hF)**2 + Lambda*S))
    err[i] = snr(f0, fSob)
    
}

plot(Lambda_list, err, type="l", col=4, ylab="SNR", xlab="lambda")

i = which.max(err)
Lambda = Lambda_list[i]

fSob = Re(fft(yF * hF / (abs(hF)**2 + Lambda*S), inverse=TRUE)) / length(yF * hF / (abs(hF)**2 + Lambda*S))

imageplot(fSob, paste('Sobolev deconvolution, SNR = ', round(snr(f0, fSob),2), 'dB'))