options(repr.plot.width=6, repr.plot.height=6)
matplot(lambda_list / lmax, W, type="l", xlab="lambda / |X * y|_inf", lty=1)
lines(c(lambda_list[i] / lmax, lambda_list[i] / lmax), c(min(W), max(W)), type="l", lty=2, col="red")
legend("right", legend=class_names[1:(length(class_names) - 3)], col = 1:(length(class_names) - 3), pch="-")