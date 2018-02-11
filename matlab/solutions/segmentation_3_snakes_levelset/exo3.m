if do_fast
phi = phi0; % initialization
clf;
k = 0;
for i=1:niter
    g0 = grad(phi,options);
    d = max(eps, sqrt(sum(g0.^2,3)) );
    g = g0 ./ repmat( d, [1 1 2] );
    K = d .* div( g,options );
    phi = phi + tau*K;
    % redistance the function from time to time
    if mod(i,30)==0
        % phi = perform_redistancing(phi);
    end
    if mod(i, floor(niter/4))==0
        k = k+1;
        subplot(2,2,k);
        plot_levelset(phi);
    end
end
end
