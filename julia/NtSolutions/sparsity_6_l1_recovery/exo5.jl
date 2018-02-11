N = 600
P = Int(N/2)
Phi = PhiRand(P, N)
klist = Array{Int64,1}(round(linspace(1, P/15., 10)))
ntrials = 2000
rip_val = zeros(length(klist), 2)

for i in 1:length(klist)
    rip_val[i, 1:2] = 0
    k = klist[i]
    for j in 1:ntrials
        I = randperm(N)
        I = I[1:k]
        (a, b) = ric(Phi[:,I])
        rip_val[i, 1:2] = max(rip_val[i, 1:2] , [a, b])
    end
end

figure(figsize = (8, 5))

plot(klist, rip_val, linewidth = 2)
plot(klist, klist.*0 + sqrt(2) - 1, "r--", linewidth = 2)
legend([L"$\delta^2_k$", L"$\delta^2_k$", "0.41"])
xlabel("k")
title(@sprintf("N = %d, P = %d", N, P))

show()
