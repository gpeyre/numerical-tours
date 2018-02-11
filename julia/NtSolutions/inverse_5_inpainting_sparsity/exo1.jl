fSpars = y
energy = []
niter = 1000
for i in 1 : niter
    fSpars = SoftThreshPsi(ProjC(fSpars, Omega), lambd)
    # record the energy
    fW = PsiS(fSpars)
    energy = [energy; 1/2*vecnorm(y-Phi(fSpars, Omega))^2 + lambd*sum(abs(fW))]
end
    
figure(figsize = (7, 5))
plot(energy, linewidth = 2)
xlabel("Iteration")
ylabel("E")
show()