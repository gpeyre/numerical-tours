plt.figure(figsize = (10,10))
U = ProjI(V)
g = np.copy(f)
k=0

tau = .03
niter = 30

for i in range(1,niter+1):
    # advect
    g = W(g,tau*V)
    U = Wt(U,tau*V)

    #Z = U
    #U[:,:,0] = W(Z[:,:,0],tau*Z)
    #U[:,:,1] = W(Z[:,:,1],tau*Z)

    # diffuse
    U = U + tau*nu*DeltaV(U)
    g = g + tau*mu*Delta(g)
    # project
    U = ProjI(U)
    # additional constraints

    #display
    if i%(niter//4) == 0:
        k +=1
        imageplot(g, "Time = %f" %(i*tau), [2,2,k])
