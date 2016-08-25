klist = [10, 30, 50]
P = 200
ntrials = 200
tmin = 0
tmax = 2.5
q = 50
t = np.linspace(tmin, tmax, q)
t1 = np.linspace(tmin, tmax, 1000)
dt = (tmax-tmin)/q

for j in range(len(klist)):
    k = klist[j]
    
    # simulation    
    v = []
    for i in range(ntrials):
        v = v + [linalg.svd(random.randn(P, k)/np.sqrt(P))[1]**2]
    
    plt.figure(figsize = (10,10))
    plt.subplot(len(klist),1,j+1)
    h = np.histogram(v, t)[0]
    h = h/np.sum(h)/dt
    plt.bar(t[:-1], h, width = 1/20, color = "darkblue", edgecolor="black")
    
    # theoritical law
    beta = k/P
    a = (1 - np.sqrt(beta))**2
    b = (1 + np.sqrt(beta))**2
    z = np.sqrt(np.maximum(t1-a, np.zeros(len(t1)))*np.maximum(b-t1, np.zeros(len(t1))))/(2*np.pi*beta*t1)
    
    plt.plot(t1, z, 'r', linewidth=3)
    plt.xlim(tmin, tmax)
    plt.ylim(0, max(h)*1.05)
    plt.title("P = %d, k = %d" %(P,k))
    
    plt.show()