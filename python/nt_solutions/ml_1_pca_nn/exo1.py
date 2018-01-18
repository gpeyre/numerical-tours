D = pairwise_distances(X0,X1)
I = np.argsort(D,axis=0)
ys = y[I]
Rmax = 70
S = np.zeros( (Rmax,1) )
for R in np.arange(1,Rmax+1):
    if R==1:
        C = ys[0,:]
    else:
        h = myhist(ys[0:R,:])
        C = np.argmax(h,axis=0)+1
        # [~,C] = max(h)
    # correct classification
    S[R-1] = np.sum(C.flatten()==y1.flatten())/n1
plt.clf
plt.bar( np.arange(1,Rmax+1), S )
plt.axis('tight')
plt.xlabel('$R$')
plt.ylabel('$S$')
