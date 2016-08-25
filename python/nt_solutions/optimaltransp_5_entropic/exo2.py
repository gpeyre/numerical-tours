plt.figure(figsize = (10,10))
glist = [.1,.01,.005,.001]
niter = 300

for k in range(len(glist)):
    gamma = glist[k]
    xi = np.exp(-C/gamma)
    b = np.ones(N[1])
    
    for i in range(niter):
        a = p/(np.dot(xi,b))
        b = q /(np.dot(np.transpose(xi),a)) 
        
    Pi = np.dot(np.dot(np.diag(a),xi),np.diag(b))
    imageplot(clamp(Pi,0,np.min(1/np.asarray(N))*.3),"$\gamma=$ %.3f" %gamma, [2,2,k+1])