

MW <- M

for ( j in 1:round(log(n,2)) ){
  
  p <- round(n/2**(j-1))
  
  sel <- seq(1, p, 1)
  
  even <- seq(1, p, 2)
  
  odd <- seq(2, p, 2)
  
  # average/ difference along X
  a <- (MW[even,sel,sel, drop=F] + MW[odd,sel,sel, drop=F])/sqrt(2)
  b <- (MW[even,sel,sel, drop=F] - MW[odd,sel,sel, drop=F])/sqrt(2)
  c <- array(0, c(dim(a)[1]+dim(b)[1], p, p))
  c[1:dim(a)[1],,] <- a
  c[(dim(a)[1]+1):(dim(a)[1]+dim(b)[1]),,] <- b
  MW[sel, sel, sel] <- c
  
  # average/ difference along Y
  a <- (MW[sel,even,sel, drop=F] + MW[sel,odd,sel, drop=F])/sqrt(2)
  b <- (MW[sel,even,sel, drop=F] - MW[sel,odd,sel, drop=F])/sqrt(2)
  c <- array(0, c(p, dim(a)[2]+dim(b)[2], p))
  c[,1:dim(a)[2],] <- a
  c[,(dim(a)[2]+1):(dim(a)[2]+dim(b)[2]),] <- b
  MW[sel, sel, sel] <- c
  
  # average/ difference along Z
  a <- (MW[sel,sel,even, drop=F] + MW[sel,sel,odd, drop=F])/sqrt(2)
  b <- (MW[sel,sel,even, drop=F] - MW[sel,sel,odd, drop=F])/sqrt(2)
  c <- array(0, c(p, p, dim(a)[3]+dim(b)[3]))
  c[,,1:dim(a)[3]] <- a
  c[,,(dim(a)[3]+1):(dim(a)[3]+dim(b)[3])] <- b
  MW[sel, sel, sel] <- c
}



