


fftshift_3d <- function(arr){
  
  swap_1 <- function(arr){
    len <- dim(arr)[1]
    half <- round(len/2)
    a <- arr[(half+1):len,,, drop=F]
    b <- arr[1:half,,, drop=F]
    c <- array(0, dim(arr))
    c[1:dim(a)[1],,] <- a
    c[(dim(a)[1]+1):(dim(a)[1]+dim(b)[1]),,] <- b
    return(c)
  }
  
  swap_2 <- function(arr){
    len <- dim(arr)[2]
    half <- round(len/2)
    a <- arr[,(half+1):len,, drop=F]
    b <- arr[,1:half,, drop=F]
    c <- array(0, dim(arr))
    c[,1:dim(a)[2],] <- a
    c[,(dim(a)[2]+1):(dim(a)[2]+dim(b)[2]),] <- b
    return(c)
  }
  
  swap_3 <- function(arr){
    len <- dim(arr)[3]
    half <- round(len/2)
    a <- arr[,,(half+1):len, drop=F]
    b <- arr[,,1:half, drop=F]
    c <- array(0, dim(arr))
    c[,,1:dim(a)[3]] <- a
    c[,,(dim(a)[3]+1):(dim(a)[3]+dim(b)[3])] <- b
    return(c)
  }
  
  arr <- swap_3(swap_2(swap_1(arr)))
  
  return(arr)
}