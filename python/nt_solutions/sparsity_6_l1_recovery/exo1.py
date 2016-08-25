eta_list = np.linspace(.1, 1, 10)
ntrials = 20
mu_mean = []
mu_std = []
N = 500
for i in range(len(eta_list)):
    eta = eta_list[i]
    P = np.round(eta*N)
    c = []
    
    for k in range(ntrials):
        c = c + [mu(PhiRand(P, N))]
        
    mu_mean = mu_mean + [np.mean(c)]
    mu_std = mu_std + [np.std(c)]
    
k_mean = 1/2*(1 + 1/np.asarray(mu_mean))

plt.figure(figsize = (10,7))

plt.subplot(2,1,1)
plt.plot(eta_list, mu_mean, linewidth=2)
plt.xlabel("$\eta$")
plt.ylabel("$mu$")

plt.subplot(2,1,2)
plt.plot(np.log10(eta_list), np.log10(k_mean), linewidth=2)
plt.xlabel("$\log_{10}(\eta)$")
plt.ylabel("$\log_{10}(1/ 2(1 + 1/ \mu))$")

plt.show()