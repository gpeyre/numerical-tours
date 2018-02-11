plt.figure(figsize = (10,5))
phi0 = np.minimum(phi1, phi2)

plt.subplot(1,2,1)
plot_levelset(phi0)
plt.title("Union")

plt.subplot(1,2,2)
plot_levelset(np.maximum(phi1, phi2))
plt.title("Intersection")

plt.show()