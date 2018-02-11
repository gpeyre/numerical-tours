fig = figure(figsize=(15,15))
#ax = gca(projection="3d")
niter = 150
stress = []
Xstress = copy(X)
ndisp = [1 5 10 min(niter,100) Inf]
k = 1
clf()
for i in 1:niter

    if ndisp[k] == i
        ax = fig[:add_subplot](2,2,k, projection = "3d")

        #swiss roll
        scatter3D(Xstress[1,:], Xstress[2,:], Xstress[3,:], s=ms, c=get_cmap("jet")((X[1,:].^2+X[3,:].^2)/100), lw=0, alpha=1)

        #graph
        I,J = findn(A)
        xx = [Xstress[1,I] Xstress[1,J]]
        yy = [Xstress[2,I] Xstress[2,J]]
        zz = [Xstress[3,I] Xstress[3,J]]

        for i in 1:length(I)
            plot(xx[i,:], yy[i,:], zz[i,:], color="black")
        end

        #params
        axis("off")
        xlim(minimum(Xstress[1,:]),maximum(Xstress[1,:]))
        ylim(minimum(Xstress[2,:]),maximum(Xstress[2,:]))
        zlim(minimum(Xstress[3,:]),maximum(Xstress[3,:]))
        ax[:view_init](elev, azim)
        k += 1
    end
    # Compute the distance matrix.
    D1 = repeat(sum(Xstress.^2, 1), outer=(n,1))
    D1 = D1 + D1' - 2*Xstress'*Xstress
    D1[D1.<0] = NaN
    D1=sqrt(D1);

    # Compute the scaling matrix.
    B = -D./max.(D1,1e-10*ones(size(D1)))
    B = B - diagm(vec(sum(B,1)))

    # update
    Xstress = (B*Xstress')'/n
    # Xstress = Xstress-repmat(mean(Xstress,2), [1 n]);
    # record stress
    D2 = abs(vec(D)-vec(D1)).^2
    append!(stress,sqrt(sum(D2[!isnan(D2)])/n^2))
end
