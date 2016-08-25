fig = plt.figure(figsize = (15,15))
klist = [1,2,4,8]
i = 0
f1 = np.copy(f)

for k in range(1,max(klist)+1):
    f1 = tW.dot(f1)
    if k == klist[i]:
        ax = fig.add_subplot(2,2,i+1,projection='3d')
        my_cmap = (np.repeat(f1[:,np.newaxis],3,1))
        ax.scatter(X0[0,:], X0[1,:], X0[2,:], lw=0, c=my_cmap, s=30)
        ax.axis("off")
        ax.view_init(elev=90, azim=-90)
        ax.dist = 6
        i = i + 1
plt.show()