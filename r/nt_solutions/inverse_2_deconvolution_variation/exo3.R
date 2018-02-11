tau = 1.9 / ( 1 + Lambda * 8 / epsilon)
fTV = y
E = rep(0, niter)

for (i in 1:niter)
{
    # Compute the gradient of the smoothed TV functional.
    Gr = grad(fTV)
    d = sqrt(epsilon**2 + (Gr**2)[,,1] + (Gr**2)[,,2])
    G = -div(Gr / array(rep(d, 2), dim=c(dim(d),2)))
    # step
    e = Phi(fTV, h) - y
    fTV = fTV - tau * ( Phi(e, h) + Lambda * G)
    # energy
    E[i] = 1/2 * norm(e)**2 + Lambda * sum(d)
}    

# display energy
plot(E, type="l", col=4, ylab="Energy", xlab="Iteration #")

imageplot(fTV)