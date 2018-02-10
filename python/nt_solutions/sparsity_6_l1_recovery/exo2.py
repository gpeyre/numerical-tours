N = 2000
P = N - 10
Phi = PhiRand(N, P)
s = rep(0, N)
s[1:6] = 1
I = supp(s)
k = length(I)

print(paste("N =", N, ", P =", P, ", |I| =", k))
print(paste("F(s) =", round(F(Phi, s), 2)))
print(paste("ERC(I) =", round(erc(Phi, I), 2)))
print(paste("w-ERC(s) =", round(werc(Phi, I), 2)))
print(paste("Coh(|s|) =", round(Coh(Phi, k), 2)))