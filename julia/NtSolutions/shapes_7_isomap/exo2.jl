J = eye(n) - ones(n,n)/n
K = -1/2.*J*D.^2*J

val, Xstrain = eigs(K, nev=2, which=:LR)

Xstrain = Xstrain.*repeat(sqrt(val)', outer=(n,1))
Xstrain = real.(Xstrain')

#plot size
figure(figsize = (15,6))

#plot points
scatter(Xstrain[1,:], Xstrain[2,:], ms, c=get_cmap("jet")((X[1,:].^2+X[3,:].^2)/100), lw=0, alpha=1)

#plot vertices
I,J = findn(A)
xx = [Xstrain[1,I] Xstrain[1,J]]
yy = [Xstrain[2,I] Xstrain[2,J]]

for i in 1:length(I)
    plot(xx[i,:], yy[i,:], color="black")
end

#params
axis("off")
xlim(minimum(Xstrain[1,:]-1),maximum(Xstrain[1,:])+1)
ylim(minimum(Xstrain[2,:]-1),maximum(Xstrain[2,:])+1)

show()
