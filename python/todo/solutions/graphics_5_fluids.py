def exo1():
    """
    Display $\Ww_{\rho U}(f)$ for various values of $\rho$.
    """
    rho_list = [2 4 8 16]
    for i in 1: length(rho_list):
        rho = rho_list(i)
        imageplot(W(f, rho*U), ['\rho = ' num2str(rho)], 2, 2, i)


def exo2():
    """
    Define an iterative scheme via:
    $$ f^{(\ell+1)} = \Ww_{\rho U}(f^{(\ell)}).  $$
    Display the result $f^{(\ell)}$, which corresponds approximately to solving
    an advection equation at time $t=\ell \rho$.
    """
    rho = .25
    niter = 12*4
    k = 0
    f1 = f
    for i in 1: niter:
        f1 = W(f1, rho*U)
        if mod(i, niter/ 4) = =0
            k = k + 1
            imageplot(f1, strcat(['t = ' num2str(i*rho)]), 2, 2, k)


def exo3():
    """
    Compute the fluid dynamic by iterating these steps.
    nitialization
    
    """
    niter = 100
    V = normalize(ProjI(V))
    g = f
    clf; k = 0
    for i in 1: niter:
        % Advect
        g = W (g, tau*U)
        V = Wt(V, tau*U)
        % Diffuse
        V = V + tau*nu*Delta(V)
        g = g + tau*mu*Delta(g)
        % Project
        V = ProjI(V)
        % Additional constraints
        if 0
            g = perform_hist_eq(g, f)
            V = normalize(V)
        % Display
        if mod(i, niter/ 4) = =0
            k = k + 1
            imageplot(g, strcat(['Time = ' num2str(i*tau)]), 2, 2, k)


