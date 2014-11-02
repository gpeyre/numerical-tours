def exo1():
    """
    Display the evolution of the image on the mesh as the number of
    iterations increases.
    """
    klist = [1 2 4 8]
    i = 1
    f1 = f
    for k in 1: max(klist):
        f1 = tW*f1
        if k = =klist(i)
            options.face_vertex_color = f1(: )
            subplot(2, 2, i)
            plot_mesh(X0, F, options)
            lighting none
            i = i + 1
    options.face_vertex_color = []


def exo2():
    """
    Determine the optimal number of iterations to maximize the SNR.
    Record, for each number |i| of iteration, the SNR in |err(i)|.
    """
    X1 = X
    err = [pnoisy]
    for i in 1: 12:
        X1 = X1*tW'
        err(i + 1) = snr(X0, X1)
        if mod(i, 2) = =0
            subplot(2, 3, i/ 2)
            plot_mesh(X1, F, options)
            axis('tight'); shading('interp')
        if err(length(err)) >max(err(1: length(err)-1))
            Xbest = X1


def exo3():
    """
    Compute the linear heat diffusion.
    Monitor the denoising
    SNR |err(l)| between $X_t$ and $X_0$ at iteration index |l|.
    """
    Xt = X
    k = 0; sob = []; err = []
    for i in 1: niter:
        % step
        Xt = Xt - tau*Xt*tL'
        % error
        err(i) = snr(X0, Xt)
        if mod(i, floor(niter/ 4)) = =0
            k = k + 1
            subplot(2, 2, k)
            plot_mesh(Xt, F, options)
            shading('interp'); axis('tight')
            % title(strcat(['T = ' num2str(Tmax*k/ 4, 3)]))


def exo4():
    """
    Solve this problem for various $\mu$ on a 3D mesh.
    Draw the evolution of the SNR denoising error as a function of $\mu$.
    """
    ntests = 15
    muList = linspace(3, 15, ntests)/ 5
    errR = []
    for i in 1: ntests:
        mu = muList(i)
        A = speye(n, n) + mu*L
    for k in 1: 3:
            Xmu(k, : ) = perform_cg(A, X(k, : )')'
        errR(i) = snr(X0, Xmu)
    plot(muList, errR, '.-'); axis('tight')
    set_label('\mu', 'SNR')


