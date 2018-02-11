figure(figsize = (10, 10))
phi = copy(phi0)
k = 0

for i in 1 : niter
    gD = NtToolBox.Grad(phi, "sym",2)
    d = max(eps().*ones(n, n), sqrt(sum(gD.^2, 3)))
    g = gD./cat(3, d[:, :], d[:, :])
    G = d.*NtToolBox.Div(g[:, :, 1], g[:, :, 2], "sym", 2) - lambd*(f0 - c1).^2 + lambd*(f0 - c2).^2
    phi = (phi + tau.*G)[:, :]
    if i % 30 == 0
        phi = perform_redistancing(phi)
    end
    if i % Base.div(niter, 4) == 0
        k = k + 1
        subplot(2, 2, k)
        NtToolBox.plot_levelset(phi, 0, f0)
    end
end
