plt.figure(figsize = (10,10))
tlist = np.linspace(.3, .9, 4)

for i in range(len(tlist)):
    T = tlist[i]
    imageplot(clamp(ThreshBlock(f, T)), "T = %.1f" %T, [2, 2, i+1])