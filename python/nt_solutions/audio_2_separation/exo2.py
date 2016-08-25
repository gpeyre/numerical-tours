sel = random.permutation(n)

sel = sel[:npts]

plt.figure(figsize = (7,5))
plt.plot(P[sel,0], P[sel,1], ".", ms = 3)
plt.xlim(-5,5)
plt.ylim(-5,5)
plt.title('Transformed domain')
plt.show()