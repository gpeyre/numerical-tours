




phi0 <- pmin(phi1, phi2)

par(mfrow=c(1,2))

plot_levelset(phi0, title="Union", lw=2)

plot_levelset(pmax(phi1, phi2), title="Intersection", lw=2)