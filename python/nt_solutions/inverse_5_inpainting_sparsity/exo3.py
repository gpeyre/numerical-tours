niter = 1000
a = U*PsiS(fSpars)
E = []
for i in range(niter):
    fTI = Psi(a)
    d = y-Phi(fTI, Omega)
    E = E + [1/2*linalg.norm(d , 'fro')**2 + lambd*np.sum(abs(a))]
    # step 
    a = SoftThresh(a + tau*PsiS(Phi(d, Omega)), lambd*tau)

plt.figure(figsize=(7,5))    
plt.plot(E)
plt.show()