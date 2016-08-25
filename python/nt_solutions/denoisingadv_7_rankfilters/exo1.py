plt.figure(figsize = (10,7))

beta_list = np.linspace(0, 1, 6)

for i in range(len(beta_list)):
        beta_c = beta_list[i]
        imageplot(phi(f, beta_c), "Beta = %.1f" %beta_c, [2, 3, i+1])