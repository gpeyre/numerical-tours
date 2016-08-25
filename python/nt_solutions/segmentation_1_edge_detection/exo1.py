plt.figure(figsize=(10,10))
slist = [1,2,5,10]

for i in range(len(slist)):
    sigma = slist[i]
    imageplot(blur(f0, sigma), "$\sigma = $ %i" %sigma, [2,2,i+1])