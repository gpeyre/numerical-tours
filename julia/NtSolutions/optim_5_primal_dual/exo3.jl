y = Phi(f0)
ProxG = (f,tau) -> f + Phi(y - Phi(f))
niter = 600
ndisp = round(linspace(1,niter, 5))
ndisp = ndisp[2:end]
E = []
f = y
g = K(y)*0
f1 = f
q = 1
for i in 1:niter    
    # update
    fold = f
    g = ProxFS( g+sigma*K(f1), sigma)
    f = ProxG(  f-tau*KS(g), tau)
    f1 = f + theta * (f-fold)
    # monitor the decay of the energ
    append!(E,F(K(f)))
    if i==ndisp[q]
        subplot(2,2,q)
        imageplot(f)
        q+=1
    end
end