lambda_list = linspace(2*1e-3, .03, 40)
err = zeros( (len(lambda_list), 1) )
for i in arange(0,len(lambda_list)):
    Lambda = lambda_list[i]
    fL2 = real( ifft2( yF * hF / ( abs(hF)**2 + Lambda) ) )
    err[i] = snr(f0,fL2)
clf;
plot(lambda_list, err)
axis('tight')
xlabel('lambda')
ylabel('SNR')

i = argmax(err);
Lambda = lambda_list[i];

fL2 = real( ifft2( yF * hF / ( abs(hF)**2 + Lambda) ) )