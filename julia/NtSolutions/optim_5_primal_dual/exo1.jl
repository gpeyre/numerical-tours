niter = 200
E = []; C = []
for i in 1:niter    
    # update
    fold = f
    g = ProxFS( g+sigma*K(f1), sigma)
    f = ProxG(  f-tau*KS(g), tau)
    f1 = f + theta * (f-fold)
    # monitor the decay of the energy
    append!(E,F(K(f)))
    append!(C,snr(f0,f))
end

plot(E)
xlim(0,niter)
