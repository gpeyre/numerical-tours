def exo1():
    """
    Implement the Floyd algorithm to compute the full distance matrix
    |D|, where |D(i,j)| is the geodesic distance between
    """
    for i in 1: n:
        % progressbar(i, n)
        D = min(D, repmat(D(: , i), [1 n]) + repmat(D(i, : ), [n 1]))


def exo2():
    """
    Perform classical MDS to compute the 2D flattening.
    entered kernel
    iagonalization
    lot graph
    """
    J = eye(n) - ones(n)/ n
    K = -1/ 2 * J*(D.^2)*J
    opt.disp = 0
    [Xstrain, val] = eigs(K, 2, 'LR', opt)
    Xstrain = Xstrain .* repmat(sqrt(diag(val))', [n 1])
    Xstrain = Xstrain'
    clf; hold on
    scatter(Xstrain(1, : ), Xstrain(2, : ), ms, v, 'filled')
    plot_graph(A, Xstrain, options)
    colormap jet(256)
    axis('equal'); axis('off')


def exo3():
    """
    Perform stress minimization MDS using SMACOF to compute the 2D flattening.
    """
    niter = 150
    stress = []
    Xstress = X
    ndisp = [1 5 10 min(niter, 100) Inf]
    k = 1
    for i in 1: niter:
        if ndisp(k) = =i
            subplot(2, 2, k)
            hold on
            scatter3(Xstress(1, : ), Xstress(2, : ), Xstress(3, : ), ms, v, 'filled')
            plot_graph(A, Xstress, options)
            colormap jet(256)
            view(v1, v2); axis('equal'); axis('off')
            k = k + 1
        % Compute the distance matrix.
        D1 = repmat(sum(Xstress.^2, 1), n, 1)
        D1 = sqrt(D1 + D1' - 2*Xstress'*Xstress)
        % Compute the scaling matrix.
        B = -D./ max(D1, 1e-10)
        B = B - diag(sum(B))
        % update
        Xstress = (B*Xstress')' / n
        % Xstress = Xstress-repmat(mean(Xstress, 2), [1 n])
        % record stress
        stress(end + 1) = sqrt(sum(abs(D(: )-D1(: )).^2) / n^2)


def exo4():
    """
    Apply Isomap to a library of small images, for instance binary digits or
    faces with a rotating camera.
    o correction for this exercise.
    """


