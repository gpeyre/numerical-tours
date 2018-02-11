## Insert your code here.

plt.figure(figsize = (10,7))

f1 = f
for i in range(6):
    f1 = medfilt(f1)
    imageplot(f1, "Iteration %i" %i, [2, 3, i+1])