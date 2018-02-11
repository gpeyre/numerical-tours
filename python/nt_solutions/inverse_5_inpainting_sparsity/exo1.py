from numpy import linalg

fSpars = y
energy = []
niter = 1000
for i in range(niter):
    fSpars = SoftThreshPsi(ProjC(fSpars, Omega), lambd)
    # record the energy
    fW = PsiS(fSpars)
    energy = energy + [1/2*linalg.norm(y-Phi(fSpars, Omega),"fro")**2 + lambd*np.sum(abs(fW))]
    
plt.figure(figsize=(7,5))
plt.plot(energy, linewidth=2)
plt.xlabel("Iteration")
plt.ylabel("E")
plt.show()