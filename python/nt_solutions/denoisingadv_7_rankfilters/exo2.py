plt.figure(figsize = (5,5))
closingopening = lambda f: opening(closing(f))
imageplot(closingopening(f))