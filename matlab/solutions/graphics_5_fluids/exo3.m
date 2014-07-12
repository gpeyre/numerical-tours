niter = 100;
% initialization
V = normalize(ProjI(V));
g = f;
% 
clf; k = 0;
for i=1:niter
    % Advect
    g = W (g,tau*U);
    V = Wt(V,tau*U);
    % Diffuse
    V = V + tau*nu*Delta(V);
    g = g + tau*mu*Delta(g);
    % Project
    V = ProjI(V);
    % Additional constraints
    if 0
        g = perform_hist_eq(g, f);
        V = normalize(V);
    end
    % Display
    if mod(i,niter/4)==0
        k = k+1;
        imageplot(g, strcat(['Time=' num2str(i*tau)]), 2,2,k);
    end
end
