def exo1():
    """
    Compute the full set of mean coordinates.
    """
    C = zeros(n, n, k)
    for i in 1: k:
        vi = V(: , i)
        U = repmat(vi, [1 n^2])-W
        nb = normalize(U)
        % length
        d = sqrt(sum(U.^2))
        s = 1
    for j in mod([i-2, i], k) + 1:
            vj = V(: , j)
            na = normalize(repmat(vj, [1 n^2])-W)
            % sign
            si = s*sign(crossp(na, nb))
            % angle
            dp = dotp(na, nb)
            theta = si .* acos(clamp(dp, -1, 1))
            % add tangent of half angle
            C(: , : , i) = C(: , : , i) + reshape(tan(theta/ 2) ./ d, [n n])
            s = -s


def exo2():
    """
    Perform the full computation of the coordinate |C(:,:,i)| by iterating the diffusion
    and imposing the boundary value.
    """
    niter = 500
    displist = round([.1 .2 .5 1]*niter); idisp = 1
    Ci = zeros(n)
    for it in 1: niter:
        Ci = (Ci(sel1, : ) + Ci(: , sel1) + Ci(sel2, : ) + Ci(: , sel2))/ 4
        Ci(I) = u
        if it = =displist(idisp)
            c = Ci + 1e-3
            c(S = =1) = 0
            B = display_shape_function(c')
            subplot(2, 2, idisp)
            hold on
            t = linspace(0, 1, n)
            imagesc(t, t, B); axis('image'); axis('off')
            colormap jet(256)
            h = plot(V(1, [1: end 1]), V(2, [1: end 1]), 'r.-')
            set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms)
            idisp = idisp + 1
    C(: , : , i) = Ci


def exo3():
    """
    Compute the full set of coordinate functions |C|.
    """
    niter = 5000
    C = zeros(n, n, k)
    for i in 1: k:
        u = zeros(k + 1, 1)
        u(i) = 1
        if i = =1
            u(end) = 1
        u = interp1(0: k, u, abscur)
        Ci = zeros(n)
    for it in 1: niter:
            Ci = (Ci(sel1, : ) + Ci(: , sel1) + Ci(sel2, : ) + Ci(: , sel2))/ 4
            Ci(I) = u
        C(: , : , i) = Ci


def exo4():
    """
    Compute the full Green Coordinates.
    """
    C = zeros(n, n, k)
    D = zeros(n, n, k)
    for i in 1: k:
        j = mod(i, k) + 1
        vi = V(: , i)
        vj = V(: , j)
        ni = N(: , i)
        %
        a = repmat(vj - vi, [1 n^2])
        b = repmat(vi, [1 n^2]) - W
        Q = sum(a.^2); s = sum(b.^2); R = 2*sum(a.*b)
        na = sqrt(sum(a.^2))
        BA = na .* sum(b .* repmat(ni, [1 n^2]))
        SRT = sqrt(4*s.*Q - R.^2)
        L0 = log(s); L1 = log(s + Q + R)
        A0 = atan(R ./ SRT) ./ SRT
        A1 = atan((2*Q + R)./ SRT) ./ SRT
        A10 = A1 - A0
        L10 = L1 - L0
        %
        d = - na .* ((4*s-(R.^2)./ Q) .* A10 + R./ (2*Q).*L10 + L1 - 2)  / (4*pi) 
        cj = + BA .* (L10./ (2*Q) - A10 .* (R./ Q)) / (2*pi)
        ci = - BA .* (L10./ (2*Q) - A10 .* (2 + R./ Q)) / (2*pi)
        D(: , : , i) = reshape(d, n, n)
        C(: , : , i) = C(: , : , i) + reshape(ci, n, n)
        C(: , : , j) = C(: , : , j) + reshape(cj, n, n)


def exo5():
    """
    Compare the Mean value, Harmonic, and Green coordinates on serveral
    cages, including a cage enclosing a caracter with two legs. Try to move
    the legs, and compare the results.
    """


def exo6():
    """
    Extend the Harmonic and Green coordinates methods to volumetric cages and volumetric data.
    """


