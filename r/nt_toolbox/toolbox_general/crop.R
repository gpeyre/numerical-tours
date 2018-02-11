

crop <- function(x, n0){
    ####
    # crop - crop an image to reduce its size
    # Only crops square black and white images for now.
    ####

    
    n = dim(x)[1]
    mid = n / 2
    # Start and end of selection
    start_ind = as.integer(ceiling(mid - (n0 / 2) + 1))
    end_ind = as.integer(ceiling(mid + (n0 / 2)))

    return (x[c(start_ind:end_ind), c(start_ind:end_ind),,])
    
  }
