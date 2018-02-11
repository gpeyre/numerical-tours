N = 2000
P = N - 10
Phi = PhiRand(N, P)
s = zeros(N)
s[1:6] = 1
I = supp(s)
k = length(I)

println(@sprintf("N = %d, P = %d, |I| = %d", N, P, k))
println(@sprintf("F(s)     = %.2f",  F(Phi, s)))
println(@sprintf("ERC(I)   = %.2f",  erc(Phi, I)))
println(@sprintf("w-ERC(s) = %.2f", werc(Phi, I)))
println(@sprintf("Coh(|s|) = %.2f", Coh(Phi, k)))
