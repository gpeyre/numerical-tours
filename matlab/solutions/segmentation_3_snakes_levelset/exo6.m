if do_fast
phi = phi0;
k = 0;
gW = grad(W,options);
for i=1:niter
    gD = grad(phi,options);
    d = max(eps, sqrt(sum(gD.^2,3)) );
    g = gD ./ repmat( d, [1 1 2] );
    G = W .* d .* div( g,options ) + sum(gW.*gD,3);
    phi = phi + tau*G;
    if mod(i,30)==0
        phi = perform_redistancing(phi);
    end
    if mod(i, floor(niter/4))==0
        k = k+1;
        subplot(2,2,k);
        plot_levelset(phi,0,f0);
    end
end
end
