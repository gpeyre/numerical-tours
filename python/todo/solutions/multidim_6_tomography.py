def exo1():
    """
    Perform the projected TV gradient desent, and monitor the epsilon-TV
    norm during the iterations.
    isplay TV decay
    """
    niter = 400
    Mtv = M1
    tv = []
    for i in 1: niter:
        % compute the gradient
        G = grad(Mtv)
        d = sqrt(sum(G.^2, 4) + epsilon^2)
        dG = -div(G ./ repmat(d, [1 1 1 3]))
        % descent 
        Mtv = Mtv - tau*dG
        % projection
        F = fftn(Mtv)
        F(mask = =1) = y
        Mtv = real(ifftn(F))
        tv(i) = sum(d(: ))
    plot(tv); axis('tight')
    xlabel('Iteration'); ylabel('TV')


