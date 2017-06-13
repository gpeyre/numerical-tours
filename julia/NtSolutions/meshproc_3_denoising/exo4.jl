ntests = 15
muList = collect(linspace(3, 15, ntests)/5.)
errR = []

for i in 1:ntests
    mu = muList[i]
    A = speye(n) + mu*L

    for k in 1:3
        Xmu[k,:] = perform_conjugate_gradient(A, X[k,:])
    end

    append!(errR, snr(X0, Xmu))
end

figure(figsize=(10,7))
plot(errR, ".-")
xlabel("mu")
ylabel("SNR")
xlim(0,ntests)
