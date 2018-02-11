G = grad(W);
G = G[:, :, 1] + 1im*G[:, :, 2]
EvalG = gamma -> NtToolBox.bilinear_interpolate(G, imag(gamma), real(gamma))
EvalW = gamma -> NtToolBox.bilinear_interpolate(W, imag(gamma), real(gamma))
#
gamma = gamma0
displist = round(linspace(0, niter, 10))
k = 1;
clf
imageplot(transpose(f))
for i in 0 : niter - 1
    n = normal(gamma)
    g = EvalW(gamma) .* normalC(gamma) - dotp(EvalG(gamma), n) .* n
    gamma = resample( gamma + dt*g )
    if i == displist[k]  
        lw = 1;    
        if (i == 0) | (i == niter)
            lw = 4;
        end
        cplot(gamma, "r", lw);
        k = k + 1;
        axis("equal"); axis("off");
    end
end
