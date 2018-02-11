niter = 100000
nsamples = 10
err_rate = 50
ElistS = zeros(div(niter,err_rate), nsamples)
for is=1:nsamples
    w = zeros(p,1)
    for it=1:niter
        if mod(it,err_rate)==1
            ElistS[1+div((it-1),err_rate),is] = E(w,X,y)
        end
        tau = tau0 / (1+it/l0)
        i = 1+Int(floor(rand()*n)) # draw uniformly
        w = w - tau * nablaEi(w,i)
    end
end
subplot(2,1,1);
#
plot(1:err_rate:niter, ElistS, "b")
plot(1+(0:length(Elist)-1)*n, Elist, "k--")
axis([1 niter minimum(ElistS) maximum(ElistS)]'); box("on")
title("E(w_l)") # ; set(gca, 'FontSize', fs);
#
subplot(2,1,2);
u = log10(ElistS-minimum(Elist))
v = log10(Elist -minimum(Elist))
plot(1:err_rate:niter, u, "b")
plot(1+(0:length(Elist)-1)*n, v, "k--")
axis([1 niter minimum(u) maximum(u)]'); box("on");
title("log(E(w_l) - min E)")#; set(gca, 'FontSize', fs);
tight_layout()