plt.figure(figsize = (5,5))
openingclosing = lambda f: closing(opening(f))
imageplot(openingclosing(f))