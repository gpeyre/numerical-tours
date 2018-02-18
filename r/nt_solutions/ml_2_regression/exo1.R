q = 50
lmax = base::norm(X0, "2")**2
lambda_list = seq(0.001, 0.3, length=q) * lmax
E = c()
W = matrix(0, q, p)

for (i in 1:q)
{
    lambda = lambda_list[i]
    w = solve(t(X0) %*% X0 + lambda * base::diag(p)) %*% t(X0) %*% y0
    W[i,] = w
    E = c(E, sqrt(sum((X1 %*% w - y1)^2 ) / n1))
}

# Display error evolution.
plot(lambda_list /lmax, E, type="l", col=4, xlab="lambda / |X|2", ylab="E")
i = which.min(E)
wRidge = W[i,]

lines(c(lambda_list[i] / lmax, lambda_list[i] / lmax), c(min(E), max(E)), type="l", lty=2, col="red")