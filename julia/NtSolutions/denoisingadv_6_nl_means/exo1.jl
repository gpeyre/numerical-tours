figure(figsize = (10, 10))

ntests = 4
tau_list = linspace(.03, .25, ntests)
e0 = -Inf
fNL = []
for k in 1 : ntests
    tau = tau_list[k]
    f1 = NLmeans(tau)
    e = snr(f0, f1)
    if e > e0
        fNL = copy(f1)
        e0 = e
    end
    imageplot(clamP(f1), @sprintf("tau = %.3f, SNR = %.1f dB" ,tau, e), [2, 2, k])
end
