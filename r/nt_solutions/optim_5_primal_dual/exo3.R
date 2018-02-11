y = Phi(f0)
ProxG = function(f,tau){f + Phi(y - Phi(f))}

niter = 600
ndisp = round(linspace(1,niter, 4))

# Update
f = y
g = K(y) * 0
f1 = f

q = 1
options(repr.plot.width=6, repr.plot.height=6)

for (i in 1:niter)
{
    # Update
    fold = f
    g = ProxFS(g + sigma * K(f1), sigma)
    f = ProxG(f - tau * KS(g), tau)
    f1 = f + theta * (f - fold)

    if (i == ndisp[q])
    {
        imageplot(f, sbpt= c(2, 2, q))
        q = q + 1
    }
    
}
