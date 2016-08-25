plt.figure(figsize=(10,10))
slist = [4, 6, 10, 15]

for i in range(len(slist)):
    sigma = slist[i]
     
    g = grad(blur(f0, sigma))
    H = hessian(blur(f0, sigma))
    a = H[:,:,0:2] * np.repeat(g[:,:,0][:,:,np.newaxis],2,axis=2) + H[:,:,1:3] * np.repeat(g[:,:,1][:,:,np.newaxis],2,axis=2)
    
    plt.subplot(2,2,i+1)
    plot_levelset(np.sum(a*g,2), 0, f0)
    plt.title("$\sigma = $ %i" %sigma)