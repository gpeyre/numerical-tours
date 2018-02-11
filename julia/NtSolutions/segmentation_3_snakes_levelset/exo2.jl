figure(figsize = (10, 5))
phi0 = min(phi1, phi2)

subplot(1, 2, 1)
NtToolBox.plot_levelset(phi0)
title("Union")

subplot(1, 2, 2)
NtToolBox.plot_levelset(max(phi1, phi2))
title("Intersection")

show()
