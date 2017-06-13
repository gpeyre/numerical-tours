function perform_conjugate_gradient(A,y)

    """
     perform_conjugate_gradient - perform (bi)-conjugate gradient

         x = perform_conjugate_gradient(A,y);

       Solves for A*x=y.
       Works for vector x,y

       Important: the algorithm assumes that the matrix is symmetric definite
       positive.

       Copyright (c) 2007 Gabriel Peyre
    """

    niter = 100
    epsilon = 1e-5
    is_sdp = 1
    x = zeros(size(A,2),1)
    normb = epsilon
    r  = y - A*x
    p = r
    r0 = sum(vec(r).^2)


    err = [sum(r0)]
    for it in 1:niter
        # auxiliary vector
        w  = A*p

        d = sum(vec(p) .* vec(w));

        if abs(d)<eps()
            d=1
        end
        alpha = repeat( [r0 / d], outer=(size(x,1), 1) );           # found optimal alpha in line search
        x = x + alpha.*p                       # new guess
        r = r - alpha.*w                       # the residual is in fact r=b-A*x

        rprev = r0;                             # save norm of the old residual
        r0 = sum(r.^2);                         # compute norm of new residual

        append!(err, sqrt( sum(r.^2) ))

        if err[end]<normb
            break
        end

        # search direction
        beta = r0./rprev;
        p = r + repeat([beta], outer=(size(x,1), 1)).*p

    end
    return x
end
