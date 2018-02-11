q = 50
lambda_list = seq(0, 20, length=q)
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
plot(lambda_list, E, type="l", col=4, xlab="lambda", ylab="E")