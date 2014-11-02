def exo1():
    """
    Compute the fixed positions $Z$ of the points indexed by $B$ that are along a
    square. Warning: $p$ is not divisible by 4.
    """
    q = floor(p/ 4)
    t = (0: (q-1))/ q
    t1 = (0: (p-3*q-1))/ (p-3*q)
    Z = [[t, t*0 + 1, 1-t, t1*0]; ...
         [t*0, t, t*0 + 1, 1-t1]]
    clf; hold('on')
    hh = plot(Z(1, [1: p 1]), Z(2, [1: p 1]), '.-')
    if using_matlab()
        set_linewidth(hh, 3)
    hh = plot(Z(1, [1: p 1]), Z(2, [1: p 1]), 'r.')
    if using_matlab()
        set_markersize(hh, 20)
    axis off; axis square


def exo2():
    """
    Compute the parameterization $Y$ on a square.
    """
    R = zeros(2, n)
    R(: , B) = Z
    if using_matlab()
        Y = (L1 \ R')'
    else
        options.maxit = 300
        Y(1, : ) = perform_cg(L1, R(1, : )', options)'
        Y(2, : ) = perform_cg(L1, R(2, : )', options)'
    plot_mesh([Y; zeros(1, n)], F)
    shading faceted; axis tight


def exo3():
    """
    Shift the $B$ positions so that the eyes of the model are approximately
    horizontal.
    olve
    """
    delta = 6
    sel = [delta: p 1: delta-1]
    R = zeros(2, n); R(: , B(sel)) = Z
    Y = (L1 \ R')'
    plot_mesh([Y; zeros(1, n)], F)
    shading faceted; axis tight


def exo4():
    """
    Compute a geometry image from the boundary free parameterization
    and use it to map a checkboard texture.
    """
    q = 64
    M = zeros(q, q, 3)
    for i in 1: 3:
        M(: , : , i) = compute_triang_interp(F, Y, X(i, : ), q)
    [Y, X] = meshgrid(1: q, 1: q)
    T = mod(X + Y, 2)
    colormap(gray(256))
    plot_surf_texture(M, T)
    view(-40, 70); zoom(1.5)
    axis tight; axis square; axis off
    camlight


def exo5():
    """
    Compute the Laplacian $L$ of the mesh, the boundary $B$ and the modified
    Laplacian $L_1$.
    ompute the symmetric Laplacian matrix.
    
    
    """
    W = make_sparse(n, n)
    for i in 1: 3:
       i1 = mod(i-1, 3) + 1
       i2 = mod(i  , 3) + 1
       i3 = mod(i + 1, 3) + 1
       pp = X(: , F(i2, : )) - X(: , F(i1, : ))
       qq = X(: , F(i3, : )) - X(: , F(i1, : ))
       % normalize the vectors   
       pp = pp ./ repmat(sqrt(sum(pp.^2, 1)), [3 1])
       qq = qq ./ repmat(sqrt(sum(qq.^2, 1)), [3 1])
       % compute angles
       a = 1 ./ tan(acos(sum(pp.*qq, 1)))
       a = max(a, 1e-2); % avoid degeneracy
       W = W + make_sparse(F(i2, : ), F(i3, : ), a, n, n)
       W = W + make_sparse(F(i3, : ), F(i2, : ), a, n, n)
    d = full(sum(W, 1))
    D = spdiags(d(: ), 0, n, n)
    L = D - W
    options.verb = 0
    B = compute_boundary(F, options)
    L1 = L
    L1(B, : ) = 0
    for i in 1: length(B):
        L1(B(i), B(i)) = 1


def exo6():
    """
    Perform the parameterization of the mesh on a circle.
    
    """
    p = length(B)
    t = linspace(0, 2*pi(), p + 1)'; t(p) = []
    x0 = cos(t); y0 = sin(t)
    Rx = zeros(n, 1); Rx(B) = x0
    Ry = zeros(n, 1); Ry(B) = y0
    if using_matlab()
        x = L1 \ Rx
        y = L1 \ Ry
    else
        options.maxit = 300
        x = perform_cg(L1, Rx, options)
        y = perform_cg(L1, Ry, options)
    Y = [x'; y']
    options.lighting = 0
    plot_mesh([Y; zeros(1, n)], F)
    shading('faceted')


def exo7():
    """
    Display a high frequency function defined on the parameteric domain on
    the mesh. What do you observe ?
    isplay the function on the 2D parameteric domain.
    """
    f = 10
    v = cos(f*Y(1, : )*2*pi()) .* cos(f*Y(2, : )*2*pi())
    options.face_vertex_color = rescale(v(: ))
    subplot(1, 2, 1)
    colormap(jet(256))
    plot_mesh([Y; zeros(1, n)], F, options)
    colormap(jet(256))
    view(2)
    subplot(1, 2, 2)
    colormap(jet(256))
    plot_mesh(X, F, options)
    colormap(jet(256))


