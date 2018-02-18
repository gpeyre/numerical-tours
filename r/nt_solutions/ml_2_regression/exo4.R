q = 200
lambda_list = seq(0.001, 0.6, length=q) * lmax
E = c()
W = matrix(0, q, p)
w = rep(0, p)

for (i in 1:q)
{
    Lambda = lambda_list[i]
    #ISTA
    for (j in 1:niter)
    {
        w = ISTA(w, Lambda,tau)
    }
    W[i,] = c(w)
    E = c(E, norm(X1 %*% w - y1) / norm(y1))
}

# Display error evolution.
plot(lambda_list /lmax, E, type="l", col=4, xlab="lambda / |X * y|_inf", ylab="E")
i = which.min(E)
wSparse = W[i,]

lines(c(lambda_list[i] / lmax, lambda_list[i] / lmax), c(min(E), max(E)), type="l", lty=2, col="red")