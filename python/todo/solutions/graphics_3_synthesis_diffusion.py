def exo1():
    """
    Perform a synthesis by running a heat diffusion, starting with a random
    noise a time |T=0|. At each step of the diffusion, perform an histogram
    equation to keep the contrast of the texture.
    tep size
    umber of iteration
    """
    T = 20
    tau = .2
    niter = round(T/ tau)
    Mheat = perform_hist_eq(randn(n), x)
    clf; k = 0
    for i in 1: niter:
        % compute the gradient
        G = div(grad(Mheat))
        % descent 
        Mheat = perform_hist_eq(Mheat + tau*G, x)
        if mod(i, floor(niter/ 4)) = =0
            k = k + 1
            imageplot(Mheat, strcat(['T = ' num2str(T*k/ 4, 3)]), 2, 2, k)


def exo2():
    """
    Starting from an initial noise image, perform a total variation
    minimization. At each step of the descent, perform an histogram
    equalization so that the texture has a flat histogram.
    tep size
    void instabilities
    umber of iteration
    """
    T = 3
    tau = .005
    epsilon = 1e-5
    niter = round(T/ tau)
    Mtv = perform_hist_eq(randn(n), x)
    clf; k = 0
    for i in 1: niter:
        % compute the gradient
        G = grad(Mtv)
        d = max(epsilon, sqrt(sum(G.^2, 3)))
        G = div(G ./ repmat(d, [1 1 2]))
        % descent 
        Mtv = perform_hist_eq(Mtv + tau*G, M)
        if mod(i, floor(niter/ 4)) = =0
            k = k + 1
            Gr = grad(Mtv)
            tv = sum(sum(sqrt(sum(Gr.^2, 3)), 2), 1)
            imageplot(Mtv, strcat(['T = ' num2str(T*k/ 4, 3) ', TV = ' num2str(tv)]), 2, 2, k)


def exo3():
    """
    Perfrom a synthesis that mixes both TV minimization (to reduce the TV
    norm)
    and wavelet histogram equalization (to control the distribution of singularities). Stop the iterations when the
    synthesized image has the same TV norm as the original one.
    """


