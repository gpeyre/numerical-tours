figure(figsize = (10, 10))
phi = copy(phi0) #initialization
k = 0
index = 1

for i in 1 : niter
    g0 = NtToolBox.Grad(phi, "sym", 2)
    d = max(eps().*ones(n, n), sqrt(sum(g0.^2, 3)))
    g = g0./cat(3, d[:, :], d[:, :])
    K = d.*NtToolBox.Div(g[:, :, 1], g[:, :, 2], "sym", 2)
    phi = (phi + tau.*K)[:, :]
    index += 1
    if i % Int(niter/4.) == 0
        k = k + 1
        subplot(2, 2, k)
        NtToolBox.plot_levelset(phi)
    end
end
