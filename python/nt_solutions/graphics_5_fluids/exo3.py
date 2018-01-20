plt.figure(figsize = (10,10))
V = normalize(ProjI(V))
g = np.copy(f)
k=0
niter=12*4

for i in range(1,niter+1):
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
    if i%(niter//4) == 0:
        k +=1
        imageplot(g, "Time = %i" %(i*tau), [2,2,k])
