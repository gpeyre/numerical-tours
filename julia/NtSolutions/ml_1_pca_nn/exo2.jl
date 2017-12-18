# bounding boxes
B = maximum(abs(Z[:,1:2]))
q = 200
r = linspace(-B,B,q)
V,U = meshgrid(r,r)
z1 = [U[:] V[:]]
#test for different R
Rlist = [1 5 10 40]
for ir=1:length(Rlist)
    R=Rlist[ir]
    #
    D = distmat(Z[:,1:2],z1);
    Ds = sort(D,1)
    I = mapslices(sortperm, D, 1)
    ys = y[I]
    #
    if R==1
        C = ys[1,:]
    else
        h = mapslices(custom_hist, ys[1:R,:],1)
        C = mapslices(indmax, h,1)
    end
    C = reshape(C, (q,q))
    # maps class to color
    Cr = zeros(q,q,3)
    for i=1:k
        for a=1:3
            Cr[:,:,a] = Cr[:,:,a] + (C.==i)*col[a,i];
        end
    end
    # display
    subplot(2,2,ir)
    imshow(permutedims(Cr,[2, 1, 3])[end:-1:1,:,:], extent=[-B, B, -B, B])#, cmap = get_cmap("jet"));
    for i=1:k
        I = find(y.==i);
        plot(Z[I,1], Z[I,2], "o", c=col[:,i]*.85)#, markeredgecolor=col[:,i]*.5 )
    #    plot(Z(I,1), Z(I,2), '.', 'Color', col(:,i), 'MarkerSize', ms);
    end
    axis("tight"); axis("equal"); axis("off")
    title(string("R=",R));
end
