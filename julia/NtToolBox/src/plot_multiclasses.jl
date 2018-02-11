function plot_multiclasses(X,y; ms=8,fs=20, disp_legend=1, disp_dim=2, edge_factor=.25)
"""
 plot_multiclasses - display data for classification

   plot_multiclasses(X,y,options);

   Copyright (c) 2017 Gabriel Peyre
"""


Xm = X -> X - repeat(mean(X,1), outer=(size(X,1), 1))
Cov = X -> Xm(X)'*Xm(X)
col = [ [1 0 0]; [0 1 0]; [0 0 1]; [0 0 0]; [0 1 1]; [1 0 1]; [1 1 0]; 
    [1 .5 .5]; [.5 1 .5]; [.5 .5 1]  ]'

p,n = size(X);
# list of classes
CL = unique(y)
k = length(CL)

# dimensionality reduction
if n>disp_dim
    U,D,V = svd(Xm(X),thin=true)
    Z = Xm(X) * V
else
    Z = X
end

lgd = []
for i=1:min(k,size(col,2))
    I = find(y.==CL[i])
    if disp_dim==2
        plot(Z[I,1], Z[I,2], "o", c=col[:,i], ms=ms,label=CL[i])
    elseif disp_dim==3 
        plot3D(Z[I,1], Z[I,2], Z[I,3], "o", c=col[:,i], ms=ms, label=CL[i])
    else
        error("Works only in 2D and 3D.")
    end
    append!(lgd,CL[i]);
end

axis("tight"); axis("equal"); box("on");gca()[:grid](false)
if disp_dim==3
    gca()[:view_init](50, 70)
end
if disp_legend==1
    legend(loc="best");
end

end