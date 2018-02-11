tau = 1.8 / ( 8/epsilon )
tau = tau * 100
E = c()
niter = 50000
ndisp = round(linspace(1, niter, 4))
x = y

q = 1

# variable for the subplots position
position = 1
options(repr.plot.width=6, repr.plot.height=6)

for (i in 1:niter)
{
    E = c(E, J(x, epsilon))
    if (i > 2  && E[i] > E[i-1])
    {
       tau = tau * .8
    }
    x = x - tau * GradTV(x, epsilon)
    x = ProjH(x)

    if (i == ndisp[q])
    {
        imageplot(x, sbpt=c(2, 2, position))
        position = position + 1
        q = q + 1
    }
}