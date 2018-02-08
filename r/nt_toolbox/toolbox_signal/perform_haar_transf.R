


perform_haar_transf <- function(M, Jmin, dir){
  
  n <- dim(M)[1]
  
  if (dir==1){
    MW <- M
    
    for (j in (Jmin:round(log(n,2))) ){
      p <- round(n/2**(j-1))
      sel <- seq(1, p, 1)
      even <- seq(1, p, 2)
      odd <- seq(2, p, 2)
      
      # average/ difference along X
      a <- (MW[even,sel,sel, drop=F] + MW[odd,sel,sel, drop=F])/sqrt(2)
      b <- (MW[even,sel,sel, drop=F] - MW[odd,sel,sel, drop=F])/sqrt(2)
      c <- array(0, c(dim(a)[1]+dim(b)[1],p,p))
      c[1:dim(a)[1],,] <- a
      c[(dim(a)[1]+1):(dim(a)[1]+dim(b)[1]),,] <- b
      MW[sel,sel,sel] <- c
      
      # average/ difference along Y
      a <- (MW[sel,even,sel, drop=F] + MW[sel,odd,sel, drop=F])/sqrt(2)
      b <- (MW[sel,even,sel, drop=F] - MW[sel,odd,sel, drop=F])/sqrt(2)
      c <- array(0, c(p,dim(a)[2]+dim(b)[2],p))
      c[,1:dim(a)[2],] <- a
      c[,(dim(a)[2]+1):(dim(a)[2]+dim(b)[2]),] <- b
      MW[sel,sel,sel] <- c
      
      # average/ difference along Z
      a <- (MW[sel,sel,even, drop=F] + MW[sel,sel,odd, drop=F])/sqrt(2)
      b <- (MW[sel,sel,even, drop=F] - MW[sel,sel,odd, drop=F])/sqrt(2)
      c <- array(0, c(p,p,dim(a)[3]+dim(b)[3]))
      c[,,1:dim(a)[3]] <- a
      c[,,(dim(a)[3]+1):(dim(a)[3]+dim(b)[3])] <- b
      MW[sel,sel,sel] <- c
      
      Mf <- MW
    }
  }
  else{
    M1 <- M
    
    for ( j in (round(log(n, 2)):(Jmin)) ){
      p <- round(n/2**j)
      sel <- seq(1, p, 1)
      sel1 <- seq(1, 2*p, 1)
      selw <- seq(p+1, 2*p, 1)
      even <- seq(1, 2*p, 2)
      odd <- seq(2, 2*p, 2)
      
      # average/ difference along X
      A <- M1[sel, sel1, sel1, drop=F]
      D <- M1[selw, sel1, sel1, drop=F]
      M1[even, sel1, sel1] <- (A + D)/sqrt(2)
      M1[odd, sel1, sel1] <- (A - D)/sqrt(2)
      
      # average/ difference along Y
      A <- M1[sel1, sel, sel1, drop=F]
      D <- M1[sel1, selw, sel1, drop=F]
      M1[sel1, even, sel1] <- (A + D)/sqrt(2)
      M1[sel1, odd, sel1] <- (A - D)/sqrt(2)
      
      # average/ difference along Z
      A <- M1[sel1, sel1, sel, drop=F]
      D <- M1[sel1, sel1, selw, drop=F]
      M1[sel1, sel1, even] <- (A + D)/sqrt(2)
      M1[sel1, sel1, odd] <- (A - D)/sqrt(2)
      
      Mf <- M1
    }
  }
  return(Mf)
}