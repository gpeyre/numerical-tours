fSpars = y
energy = c()
niter = 1000
for (i in 1:niter)
{
    fSpars = SoftThreshPsi(ProjC(fSpars, Omega), lambd)
    # record the energy
    fW = PsiS(fSpars)
    energy = c(energy, 1/2 * base::norm(y - Phi(fSpars, Omega), "F")**2 + lambd * sum(abs(fW)))
}

plot(energy, type='l', col=4, xlab="Iteration", ylab="E")

imageplot(clamp(fSpars))