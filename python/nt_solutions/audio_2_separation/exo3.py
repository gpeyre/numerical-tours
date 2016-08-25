d = np.sum(P**2, 1)
rho = .1
v = np.sort(d)
I = np.argsort(d)[::-1]

#transformed points
I = I[np.arange(1,round(rho*len(I))+1)]
P1 = P[I,:]  
    
#compute Theta
nrow = np.shape(P1)[0]
Theta = np.zeros(nrow)
for i in range(nrow):
    Theta[i] = math.atan2(P1[i,1],P1[i,0])%np.pi

nbins = 200
hist = np.histogram(Theta,nbins)
h = hist[0]/np.sum(hist[0])
t = hist[1][:-1]

plt.figure(figsize = (7,5))
plt.bar(t, h, width = np.pi/nbins, color = "darkblue", edgecolor = "darkblue")
plt.xlim(0,np.pi)
plt.ylim(0,np.max(h))
plt.show()