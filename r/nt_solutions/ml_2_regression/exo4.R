options(repr.plot.width=6, repr.plot.height=6)
sigma_list = c(.05, .1, .5, 1, 5, 10)

par(mfrow=c(3,3))

for (i in 1:length(sigma_list))
{
    sigma = sigma_list[i]
    kappa = function(X,Z){exp( -distmat(X,Z)/(2*sigma^2))}
    # Regressor.
    h = solve(K + lambda * base::diag(n)) %*% y
    Y = function(x){kappa(x,X) %*% h}
    # Eval on the grid
    yn = Y(Xn)
    dim(yn) = c(q, q)
    image(t, t, yn, main=paste("Sigma :", sigma), xlab="", ylab="")
}