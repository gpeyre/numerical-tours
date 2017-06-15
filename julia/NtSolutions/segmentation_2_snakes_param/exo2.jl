gamma = copy(gamma0);
displist = round(linspace(0, niter, 10))
k = 1;
clf; 
imageplot(transpose(W));
track  = gamma
for i in 0 : niter
    N = normal(gamma);
    g = EvalW(gamma) .* normalC(gamma) - dotp(EvalG(gamma), N) .* N;
    gamma = resample( (gamma + dt*g)); 
    if i == displist[k]   
        lw = 1
        if (i == 0) | (i == niter)
            lw = 4
        end
        cplot(gamma, "r", lw);
        k = k + 1;
        axis("equal"); axis("off");
    end
end