D = pairwise_distances(X0,X1)
I = argsort(D,axis=0)
ys = y[I]
Rmax = 70
S = zeros( (Rmax,1) )
for R in arange(1,Rmax+1):
    if R==1:
        C = ys[0,:]
    else:
        h = myhist(ys[0:R,:])
        C = argmax(h,axis=0)+1
        # [~,C] = max(h)
    # correct classification
    S[R-1] = sum(C.flatten()==y1.flatten())/n1
clf
bar( arange(1,Rmax+1), S )
axis('tight')
xlabel('$R$')
ylabel('$S$')
