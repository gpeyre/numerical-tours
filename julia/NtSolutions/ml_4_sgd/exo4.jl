ElistG = zeros(div(niter,err_rate), nsamples)
tau = .002/n
nsamples = 10
err_rate = 50
for is=1:nsamples
    w = zeros(p,1)
    G = zeros(p,n) # keep track of gradients
    g = zeros(p,1)
    for it=1:niter
        if mod(it,err_rate)==1
            ElistG[1+div(it-1,err_rate),is] = E(w,X,y)
        end
        i = 1+Int(floor(rand()*n)) # draw uniformly
        g1 = nablaEi(w,i)
        # update grad
        g = g - G[:,i] + g1
        G[:,i] = g1
        #
        w = w - tau * g
    end
end
plot(1,Inf, "b"); plot(1,Inf, "r"); plot(1,Inf, "g")
plot(1:err_rate:niter, log10(ElistS-minimum(Elist)), "b", label="SGD")
plot(1:err_rate:niter, log10(ElistA-minimum(Elist)), "r", label="SGA")
plot(1:err_rate:niter, log10(ElistG-minimum(Elist)), "g", label="SAG")
axis("tight"); box("on")
title("log(E(w_l) - min E)") # set(gca, 'FontSize', fs)

# De-duplicate legend
handles, labels = gca()[:get_legend_handles_labels]()
labels_to_display, handles_to_display = [], []
for (l,h) in zip(labels, handles)
    if !(l in labels_to_display)
        append!(labels_to_display,[l])              
        append!(handles_to_display,[h])
    end
end
legend(handles_to_display,labels_to_display);