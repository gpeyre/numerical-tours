




load_image <- function(name, n=-1, flatten=1, resc=1, grayscale=1){
  ####
  # Load an image from a file, rescale its dynamic to [0,1],
  # turn it into a grayscale image and resize it to size n x n.
  ####
  f <- load.image(name)
  # turn into normalized grayscale image
  if (grayscale == 1){
    if ( (flatten==1) & (length(dim(f))>2) ){
      f <- apply(f, c(1,2), sum)
      f <- as.cimg(f)
    }
  }
  if (resc==1){
    f <- rescale(f)
  }
  # change the size of the image
  # if (n == -1){n = 512}
  if (n > 0){
    
    if (dim(f)[4]>1){
      f <- resize(f, size_x = n, size_y = n, size_z = dim(f)[3], interpolation_type = 5)
      return(f)
    }
    
    else if (dim(f)[3]==1){
      f <- resize(f, size_x = n, size_y = n, interpolation_type = 5)
    }
    
    #f <- t(f[1:n_size,1:n_size])
    f <- t(as.matrix(f))
    return(as.cimg(f))
  }
  return(as.cimg(f))
  
}
