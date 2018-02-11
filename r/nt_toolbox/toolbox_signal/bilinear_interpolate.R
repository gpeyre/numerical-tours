




bilinear_interpolate <- function(im, x, y){
  
  x0 <- round(x)
  x1 <- x0 + 1
  y0 <- round(y)
  y1 <- y0 + 1
  
  x0[x0<1] <- 1 ; x0[x0>dim(im)[2]] <- dim(im)[2]
  x1[x1<1] <- 1 ; x1[x1>dim(im)[2]] <- dim(im)[2]
  y0[y0<1] <- 1 ; y0[y0>dim(im)[1]] <- dim(im)[1]
  y1[y1<1] <- 1 ; y1[y1>dim(im)[1]] <- dim(im)[1]
  
  Ia <- c() ; Ib <- c() ; Ic <- c() ; Id <- c()
  for (i in 1:length(x0)){
    Ia <- c(Ia, im[ y0[i], x0[i] ])
    Ib <- c(Ib, im[ y1[i], x0[i] ])
    Ic <- c(Ic, im[ y0[i], x1[i] ])
    Id <- c(Id, im[ y1[i], x1[i] ])
  }
  
  wa <- (x1-x) * (y1-y)
  wb <- (x1-x) * (y-y0)
  wc <- (x-x0) * (y1-y)
  wd <- (x-x0) * (y-y0)
  
  return(wa*Ia + wb*Ib + wc*Ic + wd*Id)
  
}
