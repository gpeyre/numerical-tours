ntests = 15
muList = np.linspace(3, 15, ntests)/5.
errR = []

for i in range(ntests):
    mu = muList[i]
    A = sparse.identity(n) + mu*L
    
    for k in range(3):
        Xmu[k,:] = sparse.linalg.cg(A, X[k,:])[0].transpose()
        
    errR = errR + [snr(X0, Xmu)]
    
plt.figure(figsize=(10,7))
plt.plot(errR, '.-')
plt.xlabel("$\mu$")
plt.ylabel("SNR")
plt.show()