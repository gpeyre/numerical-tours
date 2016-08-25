plt.figure(figsize = (14,7))

f1 = f
for i in range(4):
    f1 = opening(f1)
    imageplot(f1, "Iteration %i" %i, [2, 4, i+1])
              
f1 = f
for i in range(4):
    f1 = closing(f1)
    imageplot(f1, "Iteration %i" %i, [2, 4, i+5])