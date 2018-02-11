h = zeros(n,1)
Flist = []
tau = .5
niter = 2000
for i=1:niter
    h = h - tau * nablaF(h,K,y)
    append!(Flist, F(h,K,y))
end
plot(1:niter, Flist, lw=3)
axis("tight")