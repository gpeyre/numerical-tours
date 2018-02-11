def exo1():
    """
    Perform several step of synthesis.
    """
    M = perform_hist_eq(rand(n), M0)
    niter = 6
    for it in 1: niter:
        % randomize the offsets
        sel = randperm(w*w); sel = sel(1: noffs)
        OffX = dX(sel); OffY = dY(sel)
        % project
        M1 = zeros(n)
    for j in 1: noffs:
            ofx = OffX(j)
            ofy = OffY(j)
            Xs = mod(X + ofx-1, n) + 1
            Ys = mod(Y + ofy-1, n) + 1
            P = M(Xs + (Ys-1)*n)
    for i in 1: p*p:
                d = sum(sum((P0 - repmat(P(: , : , i), [1 1 q])).^2))
                [tmp, s] = min(d)
                P(: , : , i) = P0(: , : , s)
            M1(Xs + (Ys-1)*n) = M1(Xs + (Ys-1)*n) + P
        M = perform_hist_eq(M1 / noffs, M0)
        imageplot(M, '', 2, ceil(niter/ 2), it)


def exo2():
    """
    Perform more iteration, and increase the value of |q| and |noffs|.
    """


def exo3():
    """
    Explore the influence of the parameters |w| and |q| on the quality of
    the synthesis.
    """


def exo4():
    """
    Perform the synthesis using different textures.
    """


def exo5():
    """
    Extend the algorithm to handle color textures.
    """


def exo6():
    """
    Perform the inpainting by repeating several time the projection with
    different offsets. You do not needs to average the offset for
    inpainting.
    nitialization
    """
    noffs = 1
    niter = 6
    M = M0
    I = find(mask = =1)
    M(I) = rand(length(I), 1)
    for it in 1: niter:
        % randomize the offsets
        sel = randperm(w*w); sel = sel(1: noffs)
        OffX = dX(sel); OffY = dY(sel)
        % project
        M1 = zeros(n)
    for j in 1: noffs:
            ofx = OffX(j)
            ofy = OffY(j)
            Xs = mod(X + ofx-1, n) + 1
            Ys = mod(Y + ofy-1, n) + 1
            P = M(Xs + (Ys-1)*n)
    for i in 1: p*p:
                if sum(sum(Pmask(: , : , i))) >0
                    d = sum(sum((P0 - repmat(P(: , : , i), [1 1 q])).^2))
                    [tmp, s] = min(d)
                    P(: , : , i) = P0(: , : , s)
            M1(Xs + (Ys-1)*n) = M1(Xs + (Ys-1)*n) + P
        M = M1 / noffs
        M(mask = =0) = M0(mask = =0)
        imageplot(M, '', 2, ceil(niter/ 2), it)


def exo7():
    """
    Test the inpainting with larger holes and with various B&W and color
    images.
    """


