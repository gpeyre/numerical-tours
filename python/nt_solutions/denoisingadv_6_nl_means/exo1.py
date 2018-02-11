plt.figure(figsize = (10,10))

ntests = 4
tau_list = np.linspace(.03,.25, ntests)
e0 = -float("inf")
for k in range(ntests):
    tau = tau_list[k]
    f1 = NLmeans(tau)
    e = snr(f0,f1)
    if e > e0:
        fNL = f1
        e0 = e
    imageplot(clamp(f1),"tau = %.3f, SNR = %.1f dB" %(tau,e), [2,2,k+1])