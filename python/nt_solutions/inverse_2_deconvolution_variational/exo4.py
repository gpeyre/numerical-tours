niter = 300
lambda_list = linspace(1e-6,.01,20)
tau = 1.9 / ( 1 + max(lambda_list) * 8 / epsilon)
fBest = y 
fTV = y
err = zeros((len(lambda_list),1))
for it in  arange(0,len(lambda_list)):
    # print('.',end="",flush=True)
    Lambda = lambda_list[it]
    for i in arange(0,niter):
        # Compute the gradient of the smoothed TV functional.
        Gr = grad(fTV)
        d = sqrt(epsilon**2 + sum(Gr**2, axis=2))
        G = -div(Gr / repeat3(d,2) )
        # step
        e = Phi(fTV,h)-y
        fTV = fTV - tau*( Phi(e,h) + Lambda*G)
    err[it] = snr(f0,fTV)
    if err[it]>snr(f0,fBest):
        fBest = fTV
clf
plot(lambda_list,err)
axis('tight')
xlabel('\lambda') 
ylabel('SNR')
fTV = fBest