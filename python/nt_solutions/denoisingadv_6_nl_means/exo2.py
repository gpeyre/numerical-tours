plt.figure(figsize = (10,10))

tau = .03
q_list = [10,20]
w_list = [3,6]
ind_plot = 0
for i_q in range(len(q_list)):
    for i_w in range(len(w_list)):
        
        w = w_list[i_w]
        q = q_list[i_q]
        ind_plot += 1
        
        #patch
        w1 = 2*w + 1
        
        [X,Y,dX,dY] = np.meshgrid(np.arange(1,n+1),np.arange(1,n+1),np.arange(-w,w+1),np.arange(-w,w+1))
        X = X + dX
        Y = Y + dY
        
        X[X < 1] = 2-X[X < 1] 
        Y[Y < 1] = 2-Y[Y < 1]
         
        X[X > n] = 2*n-X[X > n]
        Y[Y > n] = 2*n-Y[Y > n]
        
        I = (X-1) + (Y-1)*n
        for i in range(n//w):
            for j in range(n//w):
                I[i,j] = np.transpose(I[i,j])
                
        patch = lambda f: np.ravel(f)[I]
        P = patch(f)
        
        #PCA
        resh = lambda P0: np.transpose((np.reshape(P0, (n*n,w1*w1), order="F")))
        remove_mean = lambda Q: Q - np.tile(np.mean(Q,0),(w1*w1,1))
        
        P1 = remove_mean(resh(P))
        C = np.dot(P1,np.transpose(P1))
        [D,V] = linalg.eig(C)
        D = np.sort(D)[::-1]
        I = np.argsort(D)[::-1]
        V = V[I,:]
        
        iresh = lambda Q: np.reshape(np.transpose(Q),(n,n,d),order="F")
        descriptor = lambda f: iresh(np.dot(np.transpose(V[: ,:d]),remove_mean(resh(P))))
        
        H = descriptor(f)
        
        #NL_means
        i = [83,72]
        
        def distance_0(i,sel): 
            H1 = (H[sel[0],:,:])
            H2 = (H1[:,sel[1],:])
            return np.sum((H2 - np.tile(H[i[0],i[1],:],(len(sel[0]),len(sel[1]),1)))**2,2)/w1*w1
        
        distance = lambda i: distance_0(i, selection(i))
        kernel = lambda i,tau: normalize(np.exp(-distance(i)/(2*tau**2)))
        selection = lambda i: np.array((clamp(np.arange(i[0]-q,i[0] + q + 1), 0, n-1), clamp(np.arange(i[1]-q,i[1] + q + 1), 0, n-1)))
        
        def NLval_0(K,sel): 
            f_temp = f[sel[0],:]
            return np.sum(K*f_temp[:, sel[1]])

        NLval = lambda i, tau: NLval_0(kernel(i, tau), selection(i))
        NLmeans = lambda tau: arrayfun(lambda i1, i2: NLval([i1,i2], tau), X, Y)
        f1 = NLmeans(tau)
        
        imageplot(clamp(f1),"q = %i, w = %i, SNR = %.1f dB" %(q,w,snr(f0,f1)), [2,2,ind_plot])