def exo1():
    """
    Implement the DR iterative algorithm on |niter| iterations.
    Keep track of the evolution of the energy $J$.
    """
    energy = []
    constr = []
    tw = w0
    for i in 1: niter:
        tw_old = tw
        w = ProxG(tw, gamma)
        rw = 2*w-tw_old
        tw = (1-mu/ 2)*tw + mu/ 2*rProxJ(rw, gamma)
        % 
        energy(i) = J(w)
        constr(i) = mynorm(A(w)-r0) / mynorm(r0)
    h = plot(min(energy, energy(1)))
    set(h, 'LineWidth', 2)
    title('J(w)')
    axis tight


