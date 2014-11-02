def exo1():
    """
    Display the effect of shearing $S_a^x$ for several values of shear $a$.
    """
    alist = linspace(0, 1, 4)
    for i in 1: length(alist):
        a = alist(i)
        imageplot(clamp(shearx(f0, a)), ['a = ' num2str(a, 2)], 2, 2, i)


def exo2():
    """
    Display the effect of rotations $R_\th$ for several values of angles $\th$.
    """
    alist = linspace(0, .99*pi/ 2, 4)
    for i in 1: length(alist):
        a = alist(i)
        imageplot(clamp(rotation(f0, a)), ['\theta = ' num2str(a, 2)], 2, 2, i)


def exo3():
    """
    Display the effect of rotations $R_\th$ for several values of angles $\th$.
    """
    alist = linspace(0, .99*pi/ 2, 4)
    for i in 1: length(alist):
        a = alist(i)
        imageplot(clamp(rotation(f0, a)), ['\theta = ' num2str(a, 2)], 2, 2, i)


def exo4():
    """
    Compute the pseudo inverse reconstruction using a conjugate gradient
    descent (function |cgs|).
    """
    cb = lambda x: reshape(Phi(PhiS(reshape(x, n, m))), n*m, 1)
    niter = 10
    [y1, FLAG, RELRES, ITER, RESVEC] = cgs(cb, y(: ), 1e-10, niter)
    y1 = reshape(y1, n, m)
    fL2 = PhiS(y1)
    imageplot(clamp(fL2))


