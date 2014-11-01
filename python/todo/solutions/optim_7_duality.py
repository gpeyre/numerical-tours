def exo1():
    """
    Perform Chambolle algorithm to solve the ROF problem.
    Monitor the primal $E$ and dual $-F$ energies.
    """
    niter = 400
    Elist = []; % primal
    Flist = []; % dual
    for i in 1: niter:
        u = proxG(u - gamma * nablaF(u), gamma)
        x = y - As(u)
        % 
        Elist(end + 1) = E(x)
        Flist(end + 1) = F(u)
    Eplot = lambda e: plot(log10(e(1: end/ 2)-min(e)))
    h = plot([Elist(7: end); -Flist(7: end)]'); axis tight
    set(h, 'LineWidth', 2)
    legend('E', '-F')


