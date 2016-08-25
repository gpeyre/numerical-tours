N = 600
P = N/2
Phi = PhiRand(P, N)
klist = np.round(np.linspace(1, P/15., 10))
ntrials = 2000
rip_val = np.zeros([len(klist),2])

for i in range(len(klist)):
    rip_val[i,:2] = 0
    k = klist[i]
    for j in range(ntrials):
        I = random.permutation(N)
        I = I[:k]
        a,b = ric(Phi[:,I])
        rip_val[i,:2] = np.maximum(rip_val[i,:2] , [a,b])
            
plt.figure(figsize = (8,5))

plt.plot(klist, rip_val, linewidth=2)
plt.plot(klist, klist*0 + np.sqrt(2)-1, 'r--', linewidth=2)
plt.legend(["$\delta^2_k$","$\delta^2_k$","0.41"])
plt.xlabel('k')
plt.title("N = %d, P = %d" %(N,P))

plt.show()