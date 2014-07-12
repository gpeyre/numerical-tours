err = [];
fHeat = y;
clf; k = 0;
for i=1:niter
    G = -div(grad(fHeat));
    fHeat = fHeat - tau*G;
    if mod(i,floor(niter/4))==0
        k = k+1;
        imageplot(clamp(fHeat), strcat(['T=' num2str(T*k/4,3)]), 2,2,k );
    end
    err(i) = snr(f0,fHeat);
    if i>1
        if err(i) > max(err(1:i-1))
            fHeat0 = fHeat;
        end
    end
end
