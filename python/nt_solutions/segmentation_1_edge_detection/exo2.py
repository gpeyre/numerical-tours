plt.figure(figsize=(10,10))
t_list = np.max(d)*np.array([1./4,1./5,1./10,1./20])

for i in range(len(t_list)):
    t = t_list[i]
    imageplot(d > t, "t = %.1f" %t , [2,2,i+1])