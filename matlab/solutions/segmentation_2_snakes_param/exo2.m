gamma = gamma0;
displist = round(linspace(1,niter,10));
k = 1;
clf; hold on;
imageplot(W);
for i=1:niter    
    N = normal(gamma);
    g = EvalW(gamma).*normalC(gamma) - dotp(EvalG(gamma), N) .* N;
    gamma = resample( gamma + dt*g );    
    if i==displist(k)       
        h = plot(imag(gamma([1:end 1])),real(gamma([1:end 1])), 'r');
        if i==1 || i==niter
            set(h, 'LineWidth', 2);
        end
        k = k+1;
        drawnow;
        axis('ij'); axis('off');
    end
end
