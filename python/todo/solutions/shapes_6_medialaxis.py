def exo1():
    """
    Compute the norm of the gradient |G| modulo |nbound|. Be careful to remove the boundary
    of the shape from this indicator. Display the thresholded gradient map.
    radient
    ompute the norm of the gadient.
    emove the boundary to the skeletton.
    """
    G = grad(Q)
    G(G <-nbound/ 2) = G(G <-nbound/ 2) + nbound
    G(G >nbound/ 2) = G(G >nbound/ 2) - nbound
    G = sqrt(sum(G.^2, 3))
    M1 = perform_convolution(M, ones(3)/ 9) >.99
    G = G.*M1
    B = G >15
    A = M
    A(B = =1) = 0
    imageplot(-A)


def exo2():
    """
    Display the Skeleton obtained for different threshold values.
    """
    thresh = [10 20 50 100]
    for i in 1: 4:
        subplot(2, 2, i)
        B = G >thresh(i)
        A = M
        A(B = =1) = 0
        imageplot(-A)


