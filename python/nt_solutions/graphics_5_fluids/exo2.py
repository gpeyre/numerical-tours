plt.figure(figsize = (10,10))
rho = .25
niter = 12*4
k=0
f1=np.copy(f)
for i in range(1,niter+1):
    f1 = W(f1, rho*U)
    if i%(niter//4) == 0:
        k +=1
        imageplot(f1, "t = %i" %(i*rho), [2,2,k])