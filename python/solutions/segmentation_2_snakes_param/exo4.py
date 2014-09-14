options.order = 2;
G = grad(W, options);
G = G(:,:,1) + 1i*G(:,:,2);
EvalG = @(gamma)interp2(1:n,1:n, G, imag(gamma), real(gamma));
EvalW = @(gamma)interp2(1:n,1:n, W, imag(gamma), real(gamma));
%
gamma = gamma0;
displist = round(linspace(1,niter,10));
k = 1;
clf; hold on;
imageplot(f);
for i=1:niter
    n = normal(gamma);
    g = EvalW(gamma).*normalC(gamma) - dotp(EvalG(gamma), n) .* n;
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
