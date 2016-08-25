g = lambda C,I: np.sum(C[:,I], 1)

plt.figure(figsize=(8,5))
dlist = np.arange(1,N/20)
criter = np.zeros([len(dlist),3])

for i in range(len(dlist)):
    s = twosparse(dlist[i])
    I = (supp(s))
    criter[i,:] = [F(Phi, s), erc(Phi,I), werc(Phi,I)]
    
criter[criter < 0] = np.float("inf")

plt.plot(dlist, criter, linewidth=2)
plt.plot(dlist, dlist*0 + 1, 'k--', linewidth=2)
plt.xlim(1,np.max(dlist))
plt.ylim(np.min(criter), np.max(criter))
plt.xlabel("d")
plt.legend(["F", "ERC", "w-ERC"])
           
plt.show()