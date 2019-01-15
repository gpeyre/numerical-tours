plt.figure(figsize = (10,10))
glist = [.1,.01,.005,.001]
niter = 300

clamp = lambda x,a,b: min(max(x,a),b)

for k in range(len(glist)):
    epsilon = glist[k]
    K = np.exp(-C/epsilon)
    v = np.ones(N[1])

    for i in range(niter):
        u = a / (np.dot(K,v))
        v = b /(np.dot(np.transpose(K),u))

    P = np.dot(np.dot(np.diag(u),K),np.diag(v))
    #imageplot(clamp(Pi,0,np.min(1/np.asarray(N))*.3),"$\gamma=$ %.3f" %gamma, [2,2,k+1])
    plt.subplot(2,2,k+1)
    plt.imshow(np.clip(P,0,np.min(1/np.asarray(N))*.3));
    #"$\gamma=$ %.3f" %gamma, [2,2,k+1])
