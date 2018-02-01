



perform_wavortho_transf <- function(f, Jmin, dir, h){
  ####
  # perform_wavortho_transf - compute orthogonal wavelet transform
  
  # fw = perform_wavortho_transf(f,Jmin,dir,options);
  
  # You can give the filter in options.h.
  
  # Works in 2D only.
  
  # Copyright (c) 2014 Gabriel Peyre
  ####
  
  
  n <- dim(f)[2]
  Jmax <- as.integer(log(n, 2) - 1)
  # compute g filter
  u <- (- rep(1, length(h) - 1)) ** (1:(length(h)-1))
  # alternate +1/-1
  g <- c(0, rev(h[2:length(h)]) * u)
  
  if (dir == 1){
    ### FORWARD ###
    fW <- as.matrix(f)
    for (j in Jmax:Jmin){
      A <- fW[1:2**(j+1), 1:2**(j+1)]
      for (d in (1:2)){
        Coarse <- subsampling(cconv(A, h, d), d)
        Detail <- subsampling(cconv(A, g, d), d)
        if (d==1){
          A <- rbind(Coarse, Detail)
        }
        else{
          A <- cbind(Coarse, Detail)
        }
      }
      fW[1:2 ** (j + 1), 1:2 ** (j + 1)] <- A
    }
    return(fW)
  }
  else{
    ### BACKWARD ###
    fW <- f
    f1 <- fW
    for (j in Jmin:Jmax){
      A <- f1[1:2 ** (j + 1), 1:2 ** (j + 1)] ###
      for (d in (1:2)){
        if (d==1){
          Coarse <- A[1:2**j, , drop=F]
          Detail <- A[(2**j + 1): 2**(j + 1), , drop=F]
        }
        else{
          Coarse <- A[ , 1:2 ** j, drop=F]
          Detail <- A[ , (2 ** j + 1):2 ** (j + 1), drop=F]
        }
        Coarse <- cconv(upsampling(Coarse, d), reverse(h), d)
        Detail <- cconv(upsampling(Detail, d), reverse(g), d)
        A <- Coarse + Detail
      }
      f1[1:2 ** (j + 1), 1:2 ** (j + 1)] <- A
    }
    return(f1)
  }
}

