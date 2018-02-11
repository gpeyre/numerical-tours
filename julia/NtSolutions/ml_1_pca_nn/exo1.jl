D = distmat(X0,X1)
Ds = sort(D,1)
I = mapslices(sortperm, D, 1)
ys = y[I]
Rmax = 50
S = []
for R=1:Rmax
    if R==1
        C = ys[1,:]
    else
        h = mapslices(custom_hist, ys[1:R,:],1)
        C = mapslices(indmax, h,1)
    end
    # correct classification 
    append!(S,sum(C[:].==y1)/n1)
end
# plot(1:Rmax, S, '.-', 'MarkerSize', ms);
bar(1:Rmax, S)
axis("tight")
axis([1 Rmax minimum(S)*.99 1]')
#SetAR(1/2)
xlabel('R'); ylabel('S');
