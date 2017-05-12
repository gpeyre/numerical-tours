sel = randperm(n)

sel = sel[1:npts]

figure(figsize = (7,5))
plot(P[sel,1], P[sel,2], ".", ms = 3)
xlim(-5,5)
ylim(-5,5)
title("Transformed domain")
