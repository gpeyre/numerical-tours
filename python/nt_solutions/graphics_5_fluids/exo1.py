plt.figure(figsize = (10,10))
rho_list = [2,4,8,16]

for i in range(len(rho_list)):
    rho = rho_list[i]
    imageplot(W(f, rho*U), "rho = %i" %rho, [2,2,i+1])