





meshgrid_2d <- function(x, y) {
  
  x <- c(x); y <- c(y)
  n <- length(x)
  m <- length(y)
  
  X <- array(rep(x, each = m),  c(m, n))
  Y <- array(rep(y, times = n),  c(m, n))
  
  return(list(X = X, Y = Y))
}




meshgrid_3d <- function(x, y, z) {
  
  x <- c(x); y <- c(y); z <- c(z)
  n <- length(x)
  m <- length(y)
  l <- length(z)
  
  X <- array(rep(x, each = m, times = l),  c(m, n, l))
  Y <- array(rep(y, each = 1, times = n*l),  c(m, n, l))
  Z <- array(rep(z, each = n*m, times = 1),  c(m, n, l))
  
  return(list(X = X, Y = Y, Z = Z))
}


meshgrid_4d <- function(x, y, z, s) {
  
  x <- c(x); y <- c(y); z <- c(z); s <- c(s)
  l_x <- length(x)
  l_y <- length(y)
  l_z <- length(z)
  l_s <- length(s)
  
  X <- array(rep(rep(x, each = 1, times = l_s*l_z), each=l_y),  c(l_y, l_x, l_z, l_s))
  Y <- array(rep(y, times=l_x*l_z*l_s),  c(l_y, l_x, l_z, l_s))
  Z <- array(rep(rep(rep(z, times=l_s), each=l_x), each=l_y),  c(l_y, l_x, l_z, l_s))
  S <- array(rep(rep(s, each=l_z*l_x), each=l_y),  c(l_y, l_x, l_z, l_s))
             
  return(list(X = X, Y = Y, Z = Z, S = S))
}



meshgrid_5d <- function(x, y, z, s, u) {
  
  x <- c(x); y <- c(y); z <- c(z); s <- c(s); u <- c(u)
  l_x <- length(x)
  l_y <- length(y)
  l_z <- length(z)
  l_s <- length(s)
  l_u <- length(u)
  
  X <- array(rep(rep(x, each = 1, times = l_u*l_s*l_z), each=l_y),  c(l_y, l_x, l_z, l_s, l_u))
  Y <- array(rep(y, times=l_x*l_z*l_s*l_u),  c(l_y, l_x, l_z, l_s, l_u))
  Z <- array(rep(rep(rep(z, times=l_s*l_u), each=l_x), each=l_y),  c(l_y, l_x, l_z, l_s, l_u))
  S <- array(rep(s, each=l_z*l_y*l_x, times=l_u),  c(l_y, l_x, l_z, l_s, l_u))
  U <- array(rep(rep(u, times=1), each=l_x*l_y*l_z*l_s),  c(l_y, l_x, l_z, l_s, l_u))
  
  return(list(X = X, Y = Y, Z = Z, S = S, U = U))
}


