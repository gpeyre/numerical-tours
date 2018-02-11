gD = NtToolBox.Grad(phi, "sym", 2)
d = max(eps().*ones(n, n), sqrt(sum(gD.^2, 3)))
g = gD ./ repeat(d[:, :], inner = [1, 1, 2])
G = d.*NtToolBox.Div(g[:, :, 1], g[:, :, 2], "sym", 2) - lambd.*(f0 - c1).^2 + lambd.*(f0 - c2).^2
