# init
I = randperm(n)
I = find(y.==1)
I = I[1:k]
C = X[I,:]
niter = 16
it_dist = 0
for it=1:niter
    # NN
    D = distmat(X,C);
    yb = mapslices(indmin,D, 2);
    # display
    if it<=3 || it==niter
        it_dist = it_dist+1
        subplot(2,2,it_dist)
        for i=1:k
            I = find(yb.==i);
            plot(Z[I,1], Z[I,2], ".", c=col[:,i], ms=12)
        end
        CV = (C-repeat(mean(X,1), outer=(k,1)))*V;
        for i=1:k
            plot(CV[i,1], CV[i,2], "o", c=col[:,i]*0.9, ms=20)
        end
        axis("tight"); axis("equal"); axis("off"); title(string("Iter 3", it))
        # update centroids
        for l=1:k
            C[l,:] = mean( X[find(yb.==l),:], 1 )
        end
    end
end