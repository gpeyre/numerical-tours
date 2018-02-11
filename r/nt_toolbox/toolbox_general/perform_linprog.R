




perform_linprog <- function(A,b,c,maxit=-1,tol=1e-10){
  
  ####
  # LINPROG  This code uses the revised simplex method to solve the linear
  # programming problem: Minimize the cost c'x subject to 
  # equations Ax=b and nonnegativity x >=0:
  # 
  # f =  Min { c'x; Ax = b, x >= 0 },
  # 
  # You must define an m by n matrix A , a column vector b with 
  # m components, and a column vector c with n components.
  # 
  # The output vector x gives the minimum cost, which is the output f.
  # 
  # [x,f,itn] = linprog(A,b,c,k,maxit,tol)
  # 
  # At most "maxit" iterations (default 10*length(b)) are applied
  # and the actual number of iterations is returned in "itn".
  # 
  # If the optimal solution is unbounded or the constraints are
  # inconsistent then a diagnostic is displayed.  
  # Bland's rule is used to resolve degeneracies.  In exact
  # arithmetic cycling is not possible.  But in real life!!!
  # If x has more than 20 components it is returned in the
  # sparse format. Note that if A has many zeros it is worth
  # passing it to linprog in sparse format. 
  # 
  # Although written for teaching purposes this routine has
  # successfully solved some problems with size(A) = [50,100000]!
  # 
  # Please report any difficulties to: idc@math.canterbury.ac.nz
  # 
  # New version                  (c) I.D.Coope, 1988, 1993
  ####
  
  c <- as.vector(c)
  b <- as.vector(b)
  m <- dim(A)[1] ; n <- dim(A)[2]
  if (maxit == -1){
    maxit <- 10*m
  }
  it <- 0
  
  D <- sign(sign(b)+.5)
  D_diag <- array(0, c(length(D), length(D)))
  diag(D_diag) <- D
  D <- D_diag                         # initial (inverse) basis matrix
  A <- cbind(A, D)                    # incorporate slack/artificial variables
  B <- (n+1):(n+m)                    # initial basis
  N <- 1:n                            # non-basis
  
  phase <- 1
  xb <- abs(b)
  s <- c(rep(0, n), rep(1, m))        # supercost
  
  
  while (phase < 3){
    df <- -1
    t <- Inf
    yb <- t(D) %*% s[B]               # multipliers for Ax=b
    while (it < maxit){
      if (length(N) == 0){
        break }
                                      # no freedom for minimization
      r <- s[N] - t(A[,N])%*%yb       # reduced costs
      rmin <- min(r)                  # determine new basic variable
      q <- which.min(r)
      if (rmin >= -tol*(max(abs(s[N])) + 1)){
        break # optimal!
      }
      it <- it+1
      if (df > 0){                   # apply Bland's rule to avoid cycling
        J <- which(r < 0)
        Nq <- min(N[J])
        q <- which(N==Nq)
      }
      d <- as.vector(t((D%*%A)[,N[q]]))
      I <- which(d > tol)
      
      if (length(I) == 0){
        print("Solution is unbounded")
        it <- -it
        break
      }
      xbd <- xb[I]/d[I]
      r <- min(xbd)
      p <- which.min(xbd)
      p <- I[p]
      if (df >= 0){                   # apply Bland's rule to avoid cycling
        J <- which(xbd == r)
        Bp <- min(B[I[J]])
        p <- which(B == Bp)
      }
      xb <- xb - (r%*%d)              # CAREFUL 
      xb[p] <- r                      # update x
      df <- r*rmin                    # change in f 
      v <- (D[p,]/d[p])               
      v <- array(v, c(1,length(v)))   # row vector
      yb <- yb + (t(v)%*%(s[N[q]] - (t(d)%*%s[B])))
      d[p] <- d[p] - 1
      v <- as.vector(v)
      D <- D - (array(d,c(length(d),1)) %*% array(v, c(1, length(v))))
      t <- B[p]
      B[p] <- N[q]
      if (t >= n){
        N <- c(N[1:q], N[(q+2):length(N)])
      }
      else{
        N[q] <- t
      }
    }
    
    xb <- xb + t(D %*% (b-(A[,B]%*%t(xb)))) # iterative refinement
    I <- which(xb < 0)                # must be due to rounding error
    if (length(I) > 0){
      xb[I] <- xb[I]-xb[I]            # so correct
    }
    if ((phase == 2) || (it < 0)){
      break                           # B, xb,n,m,res=A(:,B)*xb-b
    }
    if (xb%*%s[B] > tol){
      it <- -it
      print("No feasible solution")
      break
    }
    phase <- phase+1                  # re-initialise for Phase 2
    s <- 1e+6*max(c)*s
    s[1:n] <- c
    
  }
  
  x <- rep(0, 2*n)
  x[B] <- xb
  x <- x[1:n]
  f <- t(c) %*% x
  if (it >= maxit){
    print("Too many iterations")
    it <- -it
  }
  
  return(x)
  
} 
  



  
  
  

