f0w = W(f0)
ntrials = 100
nlaunch = 20
E0 = []
E = []
dfw = lambda fw,lambd: np.sum(abs(fw) > lambd,(1,2))

for i in range(nlaunch):
    
    fw = np.tile(f0w, (ntrials,1,1)) + sigma*random.standard_normal((ntrials,n,n))
    hfw = S(fw,lambd)
    
    #quadratic error
    e = np.sum((hfw - np.tile(f0w,(ntrials,1,1)))**2,(1,2))
    E0 = np.hstack((E0,e)) 
    
    #sure error
    e = -N*sigma**2 + np.sum((hfw-fw)**2,(1,2)) + 2*sigma**2*dfw(fw,lambd)
    E = np.hstack((E,e)) 
    
v_true = np.mean(E0)
v_sure = np.mean(E)
a = v_true - 8*np.std(E0)
b = v_true + 8*np.std(E0)
t = np.linspace(a, b, 31)
mybar = lambda e: np.histogram(e[(e > a) & (e < b)], t)

plt.figure(figsize=(10,7))

plt.subplot(2,1,1)
s = mybar(E0)[0]
plt.bar(t[:-1], s, width = (b-a)/31, color = "darkblue", edgecolor="white")
plt.axvline(v_true, color = "red", linewidth = 3)

plt.subplot(2,1,2)
s = mybar(E)[0]
plt.bar(t[:-1], s, width = (b-a)/31, color = "darkblue",edgecolor="white")
plt.axvline(v_sure, color = "red", linewidth = 3)

plt.show()