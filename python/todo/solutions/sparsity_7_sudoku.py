def exo1():
    """
    Build the encoding matrix |Aenc|. Display it.
    """
    Aenc = []
    for i in 1: n:
    for j in 1: n:
            Z = zeros(n, n, n)
            Z(i, j, : ) = 1
            Aenc(end + 1, : ) = Z(: )'
    imageplot(Aenc)


def exo2():
    """
    Build the full row matrix |Arow|. Display it.
    """
    Arow = []
    for i in 1: n:
    for k in 1: n:
            Z = zeros(n, n, n)
            Z(i, : , k) = 1
            Arow(end + 1, : ) = Z(: )'
    imageplot(Arow)


def exo3():
    """
    Build the full column matrix |Acol|. Display it.
    """
    Acol = []
    for i in 1: n:
    for k in 1: n:
            Z = zeros(n, n, n)
            Z(: , i, k) = 1
            Acol(end + 1, : ) = Z(: )'
    imageplot(Acol)


def exo4():
    """
    Create the full block matrix. Display it.
    """
    Ablock = []
    for i in 1: p:
    for j in 1: p:
    for k in 1: n:
                Z = zeros(n, n, n)
                Z((1: p) + (i-1)*p, (1: p) + (j-1)*p, k) = 1
                Ablock(end + 1, : ) = Z(: )'
    imageplot(Ablock)


def exo5():
    """
    Build the full inpainting matrix |Ainp|. Display it.
    """
    Ainp = []
    for i in 1: length(I):
        Z = zeros(n, n, n)
        Z(I(i), J(i), v(i)) = 1; % double( = =k)
        Ainp(end + 1, : ) = Z(: )'
    imageplot(Ainp)


def exo6():
    """
    Implement the Soduku solver using an interger linear programming
    algorithm.
    """


def exo7():
    """
    Perform iterative projections (POCS) on the two constraints |A*Xproj(:)=1| and
    |Xproj>=0|. Display the decay of the error |norm(A*Xproj(:)-1)| in logarithmic scale.
    """
    Xproj = encode(zeros(n))
    niter = 50000
    err = []
    for i in 1: niter:
        err(end + 1) = norm(A*Xproj(: )-1, 'fro')
    	Xproj = max(projector(Xproj), 0)
    plot(log10(err/ err(1)))
    axis('tight')


def exo8():
    """
    Prove (numerically) that for this grid, the polytope of constraints
    |P={X \ A*X(:)=1 and X>=0}| is actually reduced to a singleton, which is
    the solution of the Sudoku problem.
    """


def exo9():
    """
    Try the iterative projection on convexs set (POCS) method on this grid
    (remember that you need to re-define |A| and |pA|).
    What is your conclusion ?
    ill the constraint matrix
    OCS
    heck wether this is a valid solution.
    """
    [I, J] = ind2sub([n n], find(x1(: )~ = 0)); v = x1(x1(: )~ = 0)
    Ainp = []
    for i in 1: length(I):
        Z = zeros(n, n, n)
        Z(I(i), J(i), v(i)) = 1; % double( = =k)
        Ainp(end + 1, : ) = Z(: )'
    A = [Aenc; Arow; Acol; Ablock; Ainp]
    pA = pinv(A)
    projector = lambda u: reshape(u(: ) - pA*(A*u(: )-1), [n n n])
    Xproj = encode(zeros(n))
    err = []
    for i in 1: niter:
        err(end + 1) = norm(A*Xproj(: )-1, 'fro')
    	Xproj = clamp(projector(Xproj), 0, 1)
    plot(log10(err/ err(1)))
    axis('tight')
    [tmp, xproj] = min(abs(Xproj-1), [], 3)
    Xproj1 = encode(xproj)
    disp(['Number of violated constraints: ' num2str(sum(A*Xproj1(: )~ = 1)) '.'])


def exo10():
    """
    Compute the solution using the reweighted L1 minimization.
    Track the evolution of the number of invalidated constraints as the
    algorithm iterates.
    """
    niter = 5
    u = ones(n, n, n)
    err = []
    for i in 1: niter:
        Xrw = solvel1(A*diag(u(: ))) .* u
        u = (abs(Xrw).^(1-alpha) + epsilon)
        % 
        [tmp, xrw] = min(abs(Xrw-1), [], 3)
        Xrw1 = encode(xrw)
        err(end + 1) = sum(A*Xrw1(: )~ = 1)
    [tmp, xrw] = min(abs(Xrw-1), [], 3)
    Xrw1 = encode(xrw)
    disp(['Number of violated constraints: ' num2str(sum(A*Xrw1(: )~ = 1)) '.'])


def exo11():
    """
    Try reweighted L1 on this puzzle.
    ill matrix.
    olve
    
    """
    [I, J] = ind2sub([n n], find(x1(: )~ = 0)); v = x1(x1(: )~ = 0)
    Ainp = []
    for i in 1: length(I):
        Z = zeros(n, n, n)
        Z(I(i), J(i), v(i)) = 1; % double( = =k)
        Ainp(end + 1, : ) = Z(: )'
    A = [Aenc; Arow; Acol; Ablock; Ainp]
    niter = 6
    u = ones(n, n, n)
    err = []
    for i in 1: niter:
        Xrw = solvel1(A*diag(u(: ))) .* u
        u = (abs(Xrw).^(1-alpha) + epsilon)
        % 
        [tmp, xrw] = min(abs(Xrw-1), [], 3)
        Xrw1 = encode(xrw)
        err(end + 1) = sum(A*Xrw1(: )~ = 1)
    h = plot(err, '.-')
    set(h, 'LineWidth', 2)
    set(h, 'MarkerSize', 20)
    title('Number of invalidated constraints')


def exo12():
    """
    Try other sparsity-enforcing minimization methods, such as Orthogonal
    Matching Pursuit (OMP), or iterative hard thresholding.
    """


def exo13():
    """
    Try the different methods of this tour on a large number of Sudokus.
    """


def exo14():
    """
    Try the different methods of this tour on larger Sudokus, for |n=4,5,6|.
    """


