



fftshift <- function(input_matrix, dim = -1) {
  
  rows <- dim(input_matrix)[1]    
  cols <- dim(input_matrix)[2]    
  
  swap_up_down <- function(input_matrix) {
    rows_half <- ceiling(rows/2)
    return(rbind(input_matrix[((rows_half+1):rows), (1:cols)], input_matrix[(1:rows_half), (1:cols)]))
  }
  
  swap_left_right <- function(input_matrix) {
    cols_half <- ceiling(cols/2)
    return(cbind(input_matrix[1:rows, ((cols_half+1):cols)], input_matrix[1:rows, 1:cols_half]))
  }
  
  if (dim == -1) {
    input_matrix <- swap_up_down(input_matrix)
    return(swap_left_right(input_matrix))
  }
  else if (dim == 1) {
    return(swap_up_down(input_matrix))
  }
  else if (dim == 2) {
    return(swap_left_right(input_matrix))
  }
  else {
    stop("Invalid dimension parameter")
  }
}






fftshift_4d <- function(input_matrix, dim = -1) {
  
  rows <- dim(input_matrix)[1]    
  cols <- dim(input_matrix)[2]    
  
  swap_up_down <- function(input_matrix) {
    rows_half <- ceiling(rows/2)
    l <- length((rows_half+1):rows)
    out <- array(0, dim(input_matrix))
    out[1:l,,,] <- input_matrix[((rows_half+1):rows),,,]
    out[(l+1):dim(out)[1],,,] <- input_matrix[(1:rows_half),,,]
    return(out)
  }
  
  swap_left_right <- function(input_matrix) {
    cols_half <- ceiling(cols/2)
    l <- length((cols_half+1):cols)
    out <- array(0, dim(input_matrix))
    out[,1:l,,] <- input_matrix[, ((cols_half+1):cols),,]
    out[,(l+1):dim(out)[2],,] <- input_matrix[, 1:cols_half,,]
    return(out)
  }
  
  if (dim == -1) {
    input_matrix <- swap_up_down(input_matrix)
    return(swap_left_right(input_matrix))
  }
  else if (dim == 1) {
    return(swap_up_down(input_matrix))
  }
  else if (dim == 2) {
    return(swap_left_right(input_matrix))
  }
  else {
    stop("Invalid dimension parameter")
  }
}


