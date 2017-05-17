figure(figsize = (10,10))
rho = .25
niter = 12*4
k=0
f1=copy(f)
for i in 1:niter
    f1 = W(f1, rho*U)
    if i%(niter//4) == 0
        k +=1
        irho = i*rho
        imageplot(f1, "t = $irho", [2,2,k])
    end
end
