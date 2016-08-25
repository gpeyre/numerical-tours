f = f0 + sigma*random.standard_normal((n,n))
ntrials = 30
lambda_list = np.linspace(.3, 2, ntrials)
E0 = []
E = []

for i in range(ntrials):
    lambd = lambda_list[i]
    hf = h(f, lambd)
    E0 += [linalg.norm(hf-f0, 'fro')**2]
    E  += [SURE(f, hf, lambd)]

plt.figure(figsize = (7,5))

plt.plot(lambda_list, E0)
plt.plot(lambda_list, E)
plt.legend(["$E_0$","$E$"])
plt.xlabel("$\lambda$")
plt.xlim(np.min(lambda_list), np.max(lambda_list))
plt.ylim([np.min(E),400])

plt.show()