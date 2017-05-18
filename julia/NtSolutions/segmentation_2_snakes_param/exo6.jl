G = grad(W);
G = G[:, :, 1] + 1im*G[:, :, 2]
EvalG = gamma -> NtToolBox.bilinear_interpolate(G, imag(gamma), real(gamma))
EvalW = gamma -> NtToolBox.bilinear_interpolate(W, imag(gamma), real(gamma))
#
gamma = gamma0[:];
displist = round(linspace(0, niter, 10))
k = 1;
clf
imageplot(transpose(f))
for i in 0 : niter
    N = normal(gamma)
    g = EvalW(gamma) .* normalC(gamma) - dotp(EvalG(gamma), N) .* N
    gamma = gamma + dt*g
    gamma = resample( gamma )
    # impose start/end point
    gamma[1] = x0
    gamma[end] = x1
    if i == displist[k]   
        lw = 1;    
        if (i == 0) | (i == niter)
            lw = 4;
        end
        cplot(gamma, "r", lw);
        k = k + 1;
        axis("equal"); axis("off");
        plot(real(gamma[1]), imag(gamma[1]), "b.", markersize = 20)
        plot(real(gamma[end]), imag(gamma[end]), "b.", markersize = 20);
    end
end