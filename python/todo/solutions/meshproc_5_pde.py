def exo1():
    """
    Compute the cotangent Laplacian $W$ of the mesh
    $$ W_{i,j} = \text{cot}(\al_{i,j}) + \text{cot}(\be_{i,j}) $$
    where $\al_{i,j}$ and $\be_{i,j}$
    are the two angles oposite the the edge $(i,j)$.
    """
    W = make_sparse(n, n)
    for i in 1: 3:
       i1 = mod(i-1, 3) + 1
       i2 = mod(i  , 3) + 1
       i3 = mod(i + 1, 3) + 1
       pp = X0(: , F(i2, : )) - X0(: , F(i1, : ))
       qq = X0(: , F(i3, : )) - X0(: , F(i1, : ))
       % normalize the vectors   
       pp = pp ./ repmat(sqrt(sum(pp.^2, 1)), [3 1])
       qq = qq ./ repmat(sqrt(sum(qq.^2, 1)), [3 1])
       % compute angles
       ang = acos(sum(pp.*qq, 1))
       W = W + make_sparse(F(i2, : ), F(i3, : ), 1 ./ tan(ang), n, n)
       W = W + make_sparse(F(i3, : ), F(i2, : ), 1 ./ tan(ang), n, n)


def exo2():
    """
    Compute the linear heat diffusion.
    """
    f = f0
    k = 1
    displist = round([.05 .1 .5 1]*niter)
    for i in 1: niter:
        % step
        f = f - tau*tL*f
        if displist(k) = =i
            subplot(2, 2, k)
            options.face_vertex_color = f(: )
            plot_mesh(X0, F, options)
            lighting none
            k = k + 1


def exo3():
    """
    Compute the linear heat diffusion by iterating this gradient descent.
    """
    Tmax = 100
    niter = ceil(Tmax/ tau)
    X = X0
    k = 1
    displist = round([.05 .1 .5 1]*niter)
    for i in 1: niter:
        % step
        X = X - tau*X*tL'
        if displist(k) = =i
            subplot(2, 2, k)
            plot_mesh(X, F)
            k = k + 1


def exo4():
    """
    Solve the wave equation PDE.
    """
    f1 = f0
    f = f0
    m = 1
    for i in 1: niter:
        % step in time
        [f, f1] = update(f, f1)
        if mod(i, niter/ 4) = =1
            % display
            rho = .4
            vmax = max(abs(f))
            a = f; a(a = =max(a)) = vmax
            a(a = =min(a)) = -vmax
            options.face_vertex_color = clamp(a, -rho*vmax, rho*vmax)
            subplot(2, 2, m)
            plot_mesh(X0, F, options)
            colormap(jet(256))
            m = m + 1


