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
    N = normal(gamma);
    g = EvalW(gamma).*normalC(gamma) - dotp(EvalG(gamma), N) .* N;
    gamma = resample( gamma + dt*g );    
    % impose start/end point
    gamma(1) = x0; gamma(end) = x1;
    if i==displist(k)   
        h = plot(imag(gamma([1:end])),real(gamma([1:end])), 'r');
        if i==1 || i==niter
            set(h, 'LineWidth', 2);
        end
        h = plot(imag(gamma([1 end])),real(gamma([1 end])), 'b.'); set(h, 'MarkerSize', 30);
        axis('ij'); axis('off');
        k = k+1;
        drawnow;
    end
end
