





clip <- function(arr, arr_min, arr_max){
  ####
  # Returns the clipped version of 'arr' in the interval [arr_min, arr_max]
  ####
  d <- dim(arr)
  arr <- as.vector(arr)
  arr[arr<arr_min] <- arr_min
  arr[arr>arr_max] <- arr_max
  arr <- array(arr, d)
  return(arr)
}


