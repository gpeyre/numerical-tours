def exo1():
    """
    Implement the DR iterative algorithm on |niter| iterations.
    Keep track of the evolution of the $\ell^1$ norm.
    """
    lun = []; err = []
    tx = zeros(n, 1)
    for i in 1: niter:
        tx = (1-mu/ 2)*tx + mu/ 2*rproxG(rproxF(tx, y), gamma)
        x = proxF(tx, y)
        lun(i) = norm(x, 1)
        err(i) = norm(y-A*x)
    plot(lun)
    axis tight


def exo2():
    """
    Test the recovery of a less sparse signal.
    What do you observe ?
    """
    s = 31
    sel = randperm(n)
    x0 = zeros(n, 1); x0(sel(1: s)) = 1
    y = A*x0
    tx = zeros(n, 1)
    for i in 1: niter:
        tx = (1-mu/ 2)*tx + mu/ 2*rproxG(rproxF(tx, y), gamma)
    x = proxF(tx, y)
    subplot(2, 1, 1)
    plot_sparse_diracs(x0)
    set_graphic_sizes([], 15)
    title('Original Signal')
    subplot(2, 1, 2)
    plot_sparse_diracs(x)
    set_graphic_sizes([], 15)
    title('Recovered by L1 minimization')


def exo3():
    """
    Perform DR on the set of signals |x0|. Note that the proximal mappings
    operate in parallel on all the signals in |x0|.
    Each |i|, count the average number |proba(i)|
    of recovered vectors of sparsity |slist(i)| (up to a given, small, precision).
    """
    lun = []
    tx = x0
    niter = 2000
    for i in 1: niter:
        % progressbar(i, niter)
        tx = (1-mu/ 2)*tx + mu/ 2*rproxG(rproxF(tx, y), gamma)
        x = proxF(tx, y)
        lun(i) = norm(x, 1)
    proba = []
    E = mean(abs(x-x0)) <.05
    for j in 1: length(slist):
        s = slist(j)
        proba(j) = mean(E(Slist = =s))
    h = plot(slist, proba, 'k')
    if using_matlab()
        set(h, 'LineWidth', 2)
    axis([min(slist) max(slist) -.03 1.03])


