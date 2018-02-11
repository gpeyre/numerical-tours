q = int(p/4)
t = np.arange(0,q)/q
t1 = np.arange(0,p-3*q)/(p-3*q)
Z = np.vstack((np.hstack((t, t*0 + 1, 1-t, t1*0)),np.hstack((t*0, t, t*0 + 1, 1-t1))))
ind = np.hstack((np.arange(0,p),[0]))

plt.figure(figsize=(10,10))
plt.axis("off")
plt.xlim(-.1,1.1)
plt.ylim(-.1,1.1)
plt.plot(Z[0,ind], Z[1,ind], '.-', c="blue", linewidth=2, markerfacecolor="red", markeredgecolor="red", ms=10)
plt.show()