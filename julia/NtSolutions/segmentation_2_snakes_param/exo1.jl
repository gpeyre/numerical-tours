gamma = gamma1
displist = round(linspace(0,niter,10))
k = 1
for i in 0 : niter
    gamma = resample( gamma + dt .* normalC(gamma) )
    if i == displist[k]
        lw = 1
        if (i == 0) | (i == niter)
            lw = 4
        end
        cplot(gamma, "r", lw)
        k = k + 1
        axis("tight")
        axis("off")
    end
end