err = [];
fTV = y;
clf; k = 0;
for i=1:niter
    Gr = grad(fTV);
    d = sqrt(sum3(Gr.^2,3));
    G = -div( Gr ./ repmat( sqrt( epsilon^2 + d.^2 ) , [1 1 2]) );
    fTV = fTV - tau*G;
    if mod(i,floor(niter/4))==0
        k = k+1;
        imageplot(clamp(fTV), strcat(['T=' num2str(T*k/4,3)]), 2,2,k );
    end
    err(i) = snr(f0,fTV);
    if i>1
        if err(i) > max(err(1:i-1))
            fTV0 = fTV;
        end
    end
end
