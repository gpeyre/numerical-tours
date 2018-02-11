phi0 = min(phi1,phi2);
clf;
subplot(1,2,1);
plot_levelset(min(phi1,phi2));
title('Union');
subplot(1,2,2);
plot_levelset(max(phi1,phi2));
title('Intersection');
