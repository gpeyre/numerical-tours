import numpy as np
from numpy import linalg

def perform_linprog(A,b,c,maxit=-1,tol=1e-10):
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
    
    c = np.ravel(c, order = "F")
    b = np.ravel(b, order = "F")
    m,n = np.shape(A)
    if maxit == -1:
        maxit = 10*m
    it=0
    
    D=np.sign(np.sign(b)+.5)
    D = np.diag(D)                              # initial (inverse) basis matrix
    A = np.hstack((A,D))                        # incorporate slack/artificial variables
    B = np.arange(n,n+m)                        # initial basis
    N = np.arange(0,n)                          # non-basis
    
    phase = 1
    xb = abs(b)
    s = np.hstack((np.zeros(n), np.ones(m)))    # supercost

    while phase < 3:
        df = -1
        t = float("inf")
        yb= np.dot(np.transpose(D),s[B])        # multipliers for Ax=b
        while (it < maxit):
            if len(N) == 0:
                break   
                                                # no freedom for minimization
            r = s[N] - np.dot(np.transpose(A[:,N]),yb)          # reduced costs
            rmin = np.min(r)                                    # determine new basic variable
            q = np.argmin(r)
            if rmin >= -tol*(linalg.norm(s[N],float("inf")) + 1):
                break # optimal!
            it = it+1
            if df >= 0:                         # apply Bland's rule to avoid cycling
                J = np.where(r<0)[0]
                Nq = np.min(N[J])
                q = np.where(N==Nq)[0]
            d = np.ravel(np.dot(D,A)[:,N[q]])
            I = np.where(d > tol)[0]
            
            if len(I) == 0:
                print("Solution is unbounded")
                it = -it
                break
            xbd=xb[I]/d[I]
            r = np.min(xbd)
            p = np.argmin(xbd)
            p = I[p]
            if df >= 0:                          # apply Bland's rule to avoid cycling
                J = np.where(xbd == r)[0]
                Bp = np.min(B[I[J]])
                p = np.where(B == Bp)[0] 
            xb= xb - np.dot(r,d)                 # CAREFUL 
            xb[p] = r                            # update x
            df=r*rmin;                           # change in f 
            v = (D[p,:]/d[p])                    # row vector
            yb= yb + np.dot(np.transpose(v),s[N[q]] - np.dot(np.transpose(d),s[B]))
            d[p] = d[p] - 1
            v = np.ravel(v)
            D = D - np.dot(np.reshape(d,(np.shape(d)[0],1)),np.reshape(v,(1,np.shape(v)[0])) )                 # update inverse basis matrix
            t = B[p]
            B[p] = N[q]

            if t >= n:
                N = np.hstack((N[:q],N[(q+1):]))
            else:
                N[q] = t
                
        xb = xb + np.dot(D,b-np.dot(A[:,B],xb))  # iterative refinement
        I = np.where(xb < 0)[0]                  # must be due to rounding error
        if len(I) > 0:
            xb[I]=xb[I]-xb[I]                      # so correct
        if phase == 2 or it < 0:
            break                                # B, xb,n,m,res=A(:,B)*xb-b

        if np.dot(np.transpose(xb),s[B]) > tol:
            it=-it
            print("No feasible solution")
            break
        phase=phase+1                            # re-initialise for Phase 2
        s=1e6*linalg.norm(c,float("inf"))*s
        s[:n]=c
    
    x = np.zeros(2*n)
    x[B]=xb
    x=x[:n]
    f=np.dot(np.transpose(c),x)
    if it >= maxit:
        print("Too many iterations")
        it=-it
    return x