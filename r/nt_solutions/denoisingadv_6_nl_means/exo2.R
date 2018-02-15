

tau <- 0.003
q_list <- c(10,20)
w_list <- c(3,6)
ind_plot <- 0
for (i_q in 1:length(q_list)){
  for (i_w in 1:length(w_list)){
    
    w <- w_list[i_w]
    q <- q_list[i_q]
    ind_plot <- ind_plot + 1
    
    #patch
    w1 <- 2*w + 1
    
    grid <- meshgrid_4d(1:n, 1:n, (-w):w, (-w):w)
    X <- grid$X ; Y <- grid$Y ; dX <- grid$Z ; dY <- grid$S
    X <- X + dX
    Y <- Y + dY
    
    X[X < 1] <- 2-X[X < 1] 
    Y[Y < 1] <- 2-Y[Y < 1]
    
    X[X > n] <- 2*n-X[X > n]
    Y[Y > n] <- 2*n-Y[Y > n]
    
    I <- (X-1) + (Y-1)*n
    for (i in 1:(n%/%w)){
      for (j in 1:(n%/%w)){
        I[i,j,,] <- t(I[i,j,,])
      }
    }
    
    patch <- function(f){ array(as.vector(f)[I+1], dim(I)) }
    
    P <- patch(f)
    
    #PCA
    resh <- function(P){ t(array(P, c(n*n,w1*w1))) }
    remove_mean <- function(Q){ Q - array(rep(apply(Q, 2, mean), each=(w1*w1)), c(w1*w1, n*n)) }
    
    P1 <- remove_mean(resh(P))
    C <- P1 %*% t(P1)
    
    eg <- eigen(C)
    D <- eg$values ; V <- eg$vectors
    D <- D[order(-D)]
    I <- order(-D)
    V <- V[I,]
    
    iresh <- function(Q){ array( t(Q), c(n,n,d) ) }
    descriptor <- function(f){ iresh( t(V[, 1:d]) %*% remove_mean(resh(P)) ) }
    
    H <- descriptor(f)
    
    #NL_means
    i <- c(83,72)
    
    distance_0 <- function(i, sel){
      H1 <- H[sel[1,]+1,,]
      H2 <- H1[,sel[2,]+1,]
      return(apply((H2 - array( rep(H[i[1]+1,i[2]+1,], each=(length(sel[1,])*length(sel[2,]))), dim(H2) ))**2, c(1,2), sum)/(w1*w1))
    }
    
    distance <- function(i){ distance_0(i, selection(i)) }
    kernel <- function(i, tau){ normalize(exp(-distance(i)/ (2*tau**2))) }
    
    selection <- function(i){
      a <- clamp((i[1]-q):(i[1]+q), 0, n-1)
      b <- clamp((i[2]-q):(i[2]+q), 0, n-1)
      return( t(array(c(a,b), c(length(a), 2))) )
    }
    
    
    NLval_0 <- function(K,sel){
      f_temp <- f[sel[1,]+1,,,]
      return( sum(K*f_temp[, sel[2,]+1]) )
    }
    
    NLval <- function(i, tau){ 
      sel <- selection(i)
      K <- kernel(i, tau)
      return(NLval_0(K, sel)) }
    
    
    NLmeans <- function(tau){ arrayfun(function(i1, i2){NLval(c(i1,i2), tau)}, X, Y) }
    
    f1 <- NLmeans(tau)
    
    imageplot(clamp(f1), paste("q =", q, ", w =", w, "SNR =", round(snr(f0,f1),4), "dB"), c(2,2,ind_plot))
    
    
  }
}