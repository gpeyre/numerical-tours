N = 600
P = Int(N/2)
Phi = PhiRand(P, N)
klist = Array{Int64,1}(round(linspace(1, P/7., 20)))
ntrials = 60
proba = zeros(length(klist), 3)

for i in 1:length(klist)
    proba[i, 1:3] = 0
    k = Int(klist[i])
    for j in 1:ntrials
        s = zeros(N)
        I = randperm(N)
        I = I[1:k]
        l = randn(k, 1)
        s[I] = l./abs(l)
        proba[i, 1] = proba[i, 1] + (F(Phi, s) .< 1)
        proba[i, 2] = proba[i, 2] + (erc(Phi, I) .< 1)
        proba[i, 3] = proba[i, 3] + (werc(Phi, I) .> 0).*(werc(Phi, I) .< 1)
    end
end

figure(figsize = (8, 5))
plot(klist, proba/ntrials, linewidth = 2)
xlabel("k")
legend(["F <1", "ERC <1", "w-ERC <1"])
title(@sprintf("N = %d, P = %d", N, P))
show()
