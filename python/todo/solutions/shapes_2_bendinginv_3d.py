def exo1():
    """
    Compute the geodesic distance matrix $\de$.
    It is going to take some of time.
     progressbar(i,N);
    """
    delta = zeros(N, N)
    for i in 1: N:
        [delta(: , i), S, Q] = perform_fast_marching_mesh(V, F, i)
    delta = (delta + delta)'/ 2


def exo2():
    """
    Perform the SMACOF iterative algorithm.
    Save in a variable |s(l)| the values of
    Stress$( X^{(\ell)} )$.
    """
    niter = 50
    stress = []
    Y = V
    ndisp = [1 5 10 niter Inf]
    s = []
    k = 1
    for i in 1: niter:
        if ndisp(k) = =i
            subplot(2, 2, k)
            plot_mesh(Y, F, options)
            axis('equal'); axis('off')
            k = k + 1
        Y = Y * B(D(Y))' / N
        % update
        % Y = Y-repmat(mean(Y, 2), [1 N])
        % record stress
        s(end + 1) = Stress(D(Y))
    axis('equal'); axis('off'); axis('ij')


def exo3():
    """
    Implement a surface retrival algorithm based on these bending invariants.
    o correction for this exercise.
    """


