plt.figure(figsize=(10,10))
slist = [4, 6, 10, 15]

for i in range(len(slist)):
    sigma = slist[i]
    plt.subplot(2,2,i+1)
    plot_levelset(delta(blur(f0, sigma)) , 0, f0)
    plt.title("$\sigma = $ %i" %sigma)