


M1 <- MWT


for (j in round(log(n,2)):1){
  p <- round(n/2**j)
  sel <- seq(1, p, 1)
  sel1 <- seq(1, 2*p, 1)
  selw <- seq(p+1, 2*p)
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
  
}