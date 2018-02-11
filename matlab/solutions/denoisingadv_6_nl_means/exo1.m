ntests = 4;
tau_list = linspace(.003, .025, ntests);
e0 = -Inf; 
clf;
for k=1:ntests
    tau = tau_list(k);
    f1 = NLmeans(tau);
    e = snr(f0,f1);
    if e>e0
        fNL = f1;
        e0 = e;
    end
    imageplot(clamp(f1), ['SNR=' num2str(e,4) 'dB'], 2, 2, k);    
end
