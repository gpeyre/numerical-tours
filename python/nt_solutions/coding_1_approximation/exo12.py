plt.figure(figsize = (8,5))

names = ["regular3", "phantom", "lena", "mandrill"]

for i in range(np.shape(fList)[2]):
    fW = perform_wavortho_transf(fList[: , : , i], Jmin, + 1, h)
    cR = np.sort(np.ravel(abs(fW)))[::-1]
    err = [e for e in linalg.norm(fList[:,:,i])**2 - np.cumsum(cR**2)]
    Err = err[10:n*n//10]
    plt.plot(np.log10(np.arange(10,n*n//10)),np.log10(Err/Err[0]), linewidth=2, label = names[i])
    
plt.title("$\log_{10}(\epsilon^2[M]/ ||f||^2)$")
plt.xlim(1,np.log10(n*n//10))
plt.ylim(-7,0)
plt.legend(loc = 3)
plt.show()