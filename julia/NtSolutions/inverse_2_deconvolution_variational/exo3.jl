tau = 1.9 / ( 1 + Lambda * 8 / epsilon)
fTV = y
E = zeros(niter, 1)
for i in 1 : niter
    # Compute the gradient of the smoothed TV functional.
    Gr = grad(fTV)
    d = sqrt(epsilon^2 + sum(Gr.^2, 3))
    G = -NtToolBox.div(Gr ./ repeat3(d) )
    # step
    e = Phi(fTV,h) - y
    fTV = fTV - tau.*( Phi(e,h) + Lambda.*G)
    # energy
    E[i] = 1/2*norm(e)^2 + Lambda*sum(d)
end
# display energy
clf;
plot(E)
axis("tight")
xlabel("Iteration #")
ylabel("Energy")
