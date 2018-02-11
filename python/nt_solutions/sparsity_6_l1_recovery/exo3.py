N = 600
P = N/2
Phi = PhiRand(P, N)
klist = np.round(np.linspace(1, P/7., 20))
ntrials = 60
proba = np.zeros([len(klist),3])

for i in range(len(klist)):
    proba[i,:3] = 0
    k = klist[i]
    for j in range(ntrials):
        s = np.zeros(N)
        I = random.permutation(N)
        I = I[:k]
        s[I] = np.sign(random.randn(k, 1))
        proba[i, 0] = proba[i, 0] + (F(Phi, s) < 1)
        proba[i, 1] = proba[i, 1] + (erc(Phi, I) < 1)
        proba[i, 2] = proba[i, 2] + (werc(Phi, I) > 0)*(werc(Phi, I) < 1)
        
plt.figure(figsize = (8,5))
plt.plot(klist, proba/ntrials, linewidth=2)
plt.xlabel("k")
plt.legend(['F <1', 'ERC <1', 'w-ERC <1'])
plt.title("N = %d, P = %d" %(N,P))
plt.show()