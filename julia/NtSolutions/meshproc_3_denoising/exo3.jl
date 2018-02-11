figure(figsize=(15,15))
Xt = copy(X)
k = 0
sob = []
err = []

for i in 1:niter
    # step
    Xt = Xt - tau*Xt*tL'
    # error
    append!(err, snr(X0,Xt))
    if i % floor(niter/4) == 0
        k = k + 1
        T = Tmax*k/4
        plot_mesh(Xt, F, sub=[2,2,k], titre="T = $T")
    end
end
