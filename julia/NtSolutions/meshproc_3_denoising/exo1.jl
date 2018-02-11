fig = figure(figsize = (15,15))
klist = [1,2,4,8]
i = 1
f1 = copy(f)

for k in 1:maximum(klist)
    f1 = tW*f1
    if k == klist[i]
        ax = fig[:add_subplot](2,2,i, projection = "3d")
        my_cmap = repeat(f1,outer=(1,3))
        scatter3D(X0[1,:], X0[2,:], X0[3,:], lw=0, c=my_cmap, s=30)
        axis("off")
        ax[:view_init](90,-90)
        ax[:dist] = 6;
        i+=1
    end
end
