W = zeros(p,k)
Elist = []
tau = .01
niter = 500
for i=1:niter
    W = W - tau * nablaE(W)
    append!(Elist, E(W))
end
plot(1:niter, Elist, lw=2)
axis("tight")
