N = 2000
P = N-10
Phi = PhiRand(N, P)
s = np.zeros(N)
s[:6] = 1
I = supp(s)
k = len(I)

print("N = %d, P = %d, |I| = %d" %(N,P,k))
print("F(s)     = %.2f"  %F(Phi, s))
print("ERC(I)   = %.2f"  %erc(Phi, I))
print("w-ERC(s) = %.2f" %werc(Phi, I))
print("Coh(|s|) = %.2f" %Coh(Phi, k))