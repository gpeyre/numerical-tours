function perform_linprog(A,b,c,maxit=-1,tol=1e-10)
    """
        LINPROG  This code uses the revised simplex method to solve the linear
               programming problem: Minimize the cost c'x subject to
               equations Ax=b and nonnegativity x >=0:

                        f =  Min { c'x; Ax = b, x >= 0 },

               You must define an m by n matrix A , a column vector b with
               m components, and a column vector c with n components.

               The output vector x gives the minimum cost, which is the output f.

                        [x,f,itn] = linprog(A,b,c,k,maxit,tol)

               At most "maxit" iterations (default 10*length(b)) are applied
               and the actual number of iterations is returned in "itn".

               If the optimal solution is unbounded or the constraints are
               inconsistent then a diagnostic is displayed.
               Bland's rule is used to resolve degeneracies.  In exact
               arithmetic cycling is not possible.  But in real life!!!
               If x has more than 20 components it is returned in the
               sparse format. Note that if A has many zeros it is worth
               passing it to linprog in sparse format.

               Although written for teaching purposes this routine has
               successfully solved some problems with size(A) = [50,100000]!

               Please report any difficulties to: idc@math.canterbury.ac.nz

               New version                  (c) I.D.Coope, 1988, 1993
    """

    c = vec(c)
    b = vec(b)
    m,n = size(A)
    if maxit == -1
        maxit = 10*m
    end

    it=0

    D= sign(sign(b)+.5)
    D = diagm(D)                       # initial (inverse) basis matrix
    A = [A D]                      # incorporate slack/artificial variables
    B = collect(n+1:n+m)                     # initial basis
    N = collect(1:n)                    # non-basis

    phase = 1
    xb = abs(b)
    s =[zeros(n); ones(m)]    # supercost

    while phase < 3
        df = -1
        t = Inf
        yb= D'*s[B]        # multipliers for Ax=b
        while (it < maxit)
            if isempty(N)
                break
            end                              # no freedom for minimization
            r = s[N] - A[:,N]'*yb          # reduced costs
            rmin,q = findmin(r)                                    # determine new basic variable
            if rmin >= -tol*(norm(s[N],Inf)+1)
                break # optimal!
            end
            it = it+1
            if df >= 0                         # apply Bland's rule to avoid cycling
                J = find(r.<0)
                Nq = minimum(N[J])
                q = find(N.==Nq)[1]
            end
            d = vec(D*A[:,N[q]])
            I = find(d .> tol)
            if isempty(I)
                print("Solution is unbounded")
                it = -it
                break
            end
            xbd=xb[I]./d[I]
            r,p = findmin(xbd)
            p = I[p]
            if df >= 0                         # apply Bland's rule to avoid cycling
                J = find(xbd .== r)
                Bp = minimum(B[I[J]])
                p = find(B .== Bp)[1]
            end
            xb= xb - r*d

            xb[p] = r                            # update x
            df=r*rmin;
            # change in f
            v = vec(D[p,:])./d[p] # row vector
            yb= yb + v.*(s[N[q]] - dot(d',s[B]))
            d[p] = d[p] - 1
            D = D - d*v'                 # update inverse basis matrix
            t = B[p]
            B[p] = N[q]

            if t > n
                deleteat!(N,q)
            else
                N[q] = t
            end

        end
        xb=xb+D*(b-A[:,B]*xb);  # iterative refinement
        I = find(xb .< 0)                  # must be due to rounding error
        if length(I) > 0
            xb[I]=xb[I]-xb[I]                      # so correct
        end

        if phase == 2 || it < 0
            break                                # B, xb,n,m,res=A(:,B)*xb-b
        end
        if sum(xb.*s[B]) > tol
            it=-it
            print("No feasible solution")
            break
        end
        phase=phase+1                            # re-initialise for Phase 2
        s=1e6*norm(c, Inf)*s
        s[1:n]=c
    end
    x = zeros(2*n)
    x[B]=xb
    x=x[1:n]
    f=c'*x
    if it >= maxit
        print("Too many iterations")
        it=-it
    end
    return x
end
