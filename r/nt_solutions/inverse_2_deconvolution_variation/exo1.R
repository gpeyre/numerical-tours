Lambda_list = seq(2 * 1e-3, 0.03, length=40)
nbiter = 4000

err = rep(0, length(Lambda_list))

for (i in 1:length(Lambda_list))
{
    Lambda = Lambda_list[i]
    fL2 = Re(fft(yF * hF / (abs(hF)**2 + Lambda), inverse=TRUE)) / length(yF * hF / (abs(hF)**2 + Lambda))
    err[i] = snr(f0, fL2)
}

plot(Lambda_list, err, type="l", col=4, ylab="SNR", xlab="lambda")

i = which.max(err)
Lambda = Lambda_list[i]

fL2 = Re(fft( yF * hF / ( abs(hF)**2 + Lambda), inverse=TRUE)) / length(yF * hF / (abs(hF)**2 + Lambda))

imageplot(fL2, paste('L2 deconvolution, SNR = ', round(snr(f0, fL2),2), 'dB'))