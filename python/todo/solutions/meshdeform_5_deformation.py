def exo1():
    """
    Perform a more complicated deformation of the boundary.
    eform
    isplay it.
    """
    Delta0 = zeros(3, n)
    Delta0(3, I) = vertex(1, I).^2 - vertex(2, I).^2
    vertex1 = vertex + (L1 \ Delta0')'
    plot_mesh(vertex1, faces)
    view(-150, 45)


def exo2():
    """
    Move both the inside and the boundary.
    
    aplacian
    eform
    isplay it.
    """
    J = reshape(1: p^2, p, p)
    q = round(.3*p); sel = q: p-q + 1
    I0 = unique([J(1, : ) J(end, : ) J(: , 1)' J(: , end)']')
    I1 = unique([J(q, sel) J(p-q + 1, sel) J(sel, q)' J(sel, p-q + 1)']')
    I = [I0; I1]
    Delta0 = zeros(3, n)
    Delta0(3, I1) = 1
    L1 = L; L1(I, : ) = 0; L1(I + (I-1)*n) = 1
    vertex1 = vertex + (L1 \ Delta0')'
    plot_mesh(vertex1, faces)
    view(-150, 45)
    zoom(.8)


def exo3():
    """
    Apply the mesh deformation method to a real mesh, with both large scale
    and fine scale details.
    """


def exo4():
    """
    Apply the deformation to the coarse mesh |vertex0| to obtain |vertex1|.
    *Important:* you need to compute and use the cotan Laplacian of the coarse
    mesh, not of the original mesh!
    ompute laplacian
    
    aplacian
    eform
    isplay it.
    """
    W0 = sparse(n, n)
    for i in 1: 3:
       i1 = mod(i-1, 3) + 1
       i2 = mod(i  , 3) + 1
       i3 = mod(i + 1, 3) + 1
       pp = vertex0(: , faces(i2, : )) - vertex0(: , faces(i1, : ))
       qq = vertex0(: , faces(i3, : )) - vertex0(: , faces(i1, : ))
       pp = pp ./ repmat(sqrt(sum(pp.^2, 1)), [3 1])
       qq = qq ./ repmat(sqrt(sum(qq.^2, 1)), [3 1])
       ang = acos(sum(pp.*qq, 1))
       u = cot(ang)
       u = clamp(u, 0.01, 100)
       W0 = W0 + sparse(faces(i2, : ), faces(i3, : ), u, n, n)
       W0 = W0 + sparse(faces(i3, : ), faces(i2, : ), u, n, n)
    d = full(sum(W0, 1))
    D = spdiags(d(: ), 0, n, n)
    L0 = D - W0
    Delta0 = zeros(3, n)
    Delta0(3, I1) = 1
    L01 = L0; L01(I, : ) = 0; L01(I + (I-1)*n) = 1
    vertex1 = vertex0 + (L01 \ Delta0')'
    plot_mesh(vertex1, faces)
    view(-150, 45)
    zoom(.8)


def exo5():
    """
    Add the normal contribution |d.*normal| to |vertex1|, but
    after replacing the normal of |vertex0| by the normal of |vertex1|.
    isplay it.
    """
    normal1 = compute_normal(vertex1, faces)
    vertex2 = vertex1 + d .* normal1
    plot_mesh(vertex2, faces)
    view(-150, 45)


def exo6():
    """
    Try on other surfaces. How can you compute |vertex0| for an arbitrary
    surface ?
    """


def exo7():
    """
    Compute the bi-laplacian deformation of the coarse shape
    |vertex0| by using |LL| instead of |L|.
    What do you observe ?
    eform
    isplay it.
    """
    LL1 = L0*L0
    LL1(I, : ) = 0
    LL1(I + (I-1)*n) = 1
    vertex1 = vertex0 + (LL1 \ Delta0')'
    plot_mesh(vertex1, faces)
    view(-150, 45)
    zoom(.8)


def exo8():
    """
    Compute the deformation obtained by moving from the Laplacian
    to the bi-laplacian, i.e. with |t*L+(1-t)*LL| for varying t.
    """
    t = [0 .02 .1 1]
    LL0 = L0*L0
    for i in 1: 4:
        LL1 = t(i)*L0 + (1-t(i))*LL0
        LL1(I, : ) = 0
        LL1(I + (I-1)*n) = 1
        % deform
        vertex1 = vertex0 + (LL1 \ Delta0')'
        % Display it.
        subplot(2, 2, i)
        plot_mesh(vertex1, faces)
        view(-150, 45)
        zoom(.8)


def exo9():
    """
    Apply the full model (Laplacian, bi-Laplacian and non-linear
    deformation) to the deformation of a complicated mesh.
    """


