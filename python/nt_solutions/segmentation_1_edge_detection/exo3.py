plt.figure(figsize=(10,10))
slist = [1, 2, 4, 6]

for i in range(len(slist)):
    sigma = slist[i]
    d = np.sqrt(np.sum(nabla(blur(f0, sigma))**2, 2))
    t = np.max(d)*1./5
    imageplot(d > t, "$\sigma =$ %.1f" %sigma , [2,2,i+1])