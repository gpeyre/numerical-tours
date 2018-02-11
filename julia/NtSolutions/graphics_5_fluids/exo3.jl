figure(figsize = (10,10))
V = Normalize(ProjI(V))
g = copy(f)
k=0

for i in 1:niter
    # advect
    g = W(g,tau*U)
    V = Wt(V,tau*U)
    # diffuse
    V = V + tau*nu*Delta(V)
    g = g + tau*mu*Delta(g)
    # project
    V = ProjI(V)
    # additional constraints

    #display
    if i%(Base.div(niter,4)) == 0
        k +=1; itau = i*tau
        imageplot(g, "Time = $itau", [2,2,k])
    end
end
