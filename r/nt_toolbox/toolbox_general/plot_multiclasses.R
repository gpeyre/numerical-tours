Xm = function(X){as.matrix(X - rep(colMeans(X), rep.int(nrow(X), ncol(X))))}

plot_multiclasses = function(X, y, dim)
{
    
    svd_decomp = svd(Xm(X))
    U = svd_decomp$u
    D = svd_decomp$d
    V = svd_decomp$v
    Z = Xm(X) %*% V
    classes = sort(unique(y))
    nb_classes = length(classes)

    if (dim == 2)
    {
        cols = c(2:(nb_classes + 1))
        for (i in classes)
        {
            
            I = (y==i)
            plot(Z[I,1], Z[I,2], col=cols[i], xlim=c(min(Z[,1]), max(Z[,1])),  
                 ylim=c(min(Z[,2]), max(Z[,2])), xlab="", ylab="", pch=16)
            par(new=TRUE)
        }

    
    legend("topright", legend=classes, col=cols, pch="o")
    }
    else if(dim == 3)
    {
        first = TRUE
        cols = c(2:(nb_classes + 1))
        for (i in 1:length(classes))
        {
            I = (y==classes[i])
            if (first==TRUE)
            {
                scatter3D(Z[I,1], Z[I,2], Z[I,3], col=cols[i], xlim=c(min(Z[,1]), max(Z[,1])),  
                          ylim=c(min(Z[,2]), max(Z[,2])), zlim=c(min(Z[,3]), max(Z[,3])), xlab="", ylab="",
                          zlab="", pch=16)
                first = FALSE
                par(new=TRUE)
                
            }
            else
            {
                 scatter3D(Z[I,1], Z[I,2], Z[I,3], col=cols[i], xlim=c(min(Z[,1]), max(Z[,1])),  
                           ylim=c(min(Z[,2]), max(Z[,2])), zlim=c(min(Z[,3]), max(Z[,3])), xlab="", ylab="",
                           zlab="", add=TRUE, pch=16)
            }
        }

        
        legend("topright", legend=classes, col=cols, pch="o")
    }

}
