def exo1():
    """
    Perform the projected gradient descent.
    Record in a variable |E| the evolution of the Sobolev energy $E$.
    """
    niter = 50
    E = []
    k = 1; ndisp = [1 5 10 niter]
    norm1 = lambda f: norm(f(: ))
    f = y
    for i in 1: niter:
        E(i) = norm1(grad(f))
        f = Pi(f + tau*Delta(f))
        if i = =ndisp(k)
            imageplot(f, ['iter = ' num2str(i)], 2, 2, k)
            k = k + 1


def exo2():
    """
    Perform the projected gradient descent.
    Record in a variable |J| the evolution of the TV energy $J_\epsilon$.
    """
    niter = 700
    J = []
    k = 1; ndisp = ceil(niter*[.001 1/ 10 1/ 5 1])
    f = y
    for i in 1: niter:
        J(i) = sum(sum(Amplitude(grad(f))))
        f = Pi(f - tau*G(f))
        if i = =ndisp(k)
            imageplot(f, ['iter = ' num2str(i)], 2, 2, k)
            k = k + 1


def exo3():
    """
    Perform Sobolev inpainting on this image.
    """
    Pi = lambda f: f.*(1-Gamma) + y.*Gamma
    tau = .8/ 4
    niter = 300
    E = []
    k = 1; ndisp = [1 5 10 niter]
    norm1 = lambda f: norm(f(: ))
    f = y; f(Gamma = =0) = .5
    for i in 1: niter:
        E(i) = norm1(grad(f))
        f = Pi(f + tau*Delta(f))
    imageplot(f)


def exo4():
    """
    Try other methods to solve this inpainting problem.
    You can for instance have a look on the numerical on sparsity for
    deconvolution and inpainting.
    """


