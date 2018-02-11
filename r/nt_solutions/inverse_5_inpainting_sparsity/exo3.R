niter = 1000

a = U * PsiS(fSpars)
E = c()

for (i in 1:niter)
{
    fTI = Psi(a)
    d = y - Phi(fTI, Omega)
    E = c(E,1/2 * base::norm(d, "F")**2 + lambd * sum(abs(a)))
    # step 
    a = SoftThresh(a + tau * PsiS(Phi(d, Omega)), lambd * tau)
}

plot(E, type='l', col=4, xlab="Iteration", ylab="E")

imageplot(clamp(fTI))