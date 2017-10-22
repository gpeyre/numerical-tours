


isosurface <- function(M,v,step,title=""){
  ####
  # returns the isosurface of value v of M, subsetting M with the steps argument
  ####

  sel <- seq(1, dim(M)[1], step)
  
  triangles <- computeContour3d(M[sel, sel, sel],level = v)
  
  
  # rotation matrix
  u <- (10/360)*(2*pi)
  v <- (0/360)*(2*pi)
  w <- (155/360)*(2*pi)
  
  R_mat <- matrix(c(cos(v)*cos(w),
                    cos(v)*sin(w),
                    -sin(v),
                    0,
                    sin(u)*sin(v)*cos(w) - cos(u)*sin(w),
                    cos(u)*cos(w) + sin(u)*sin(v)*sin(w),
                    sin(u)*cos(v),
                    0,
                    sin(u)*sin(w) + cos(u)*sin(v)*cos(w),
                    cos(u)*sin(v)*sin(w) - sin(u)*cos(w),
                    cos(u)*cos(v),
                    0,
                    0,
                    0,
                    0,
                    1), c(4,4))
  
  color_vec <- rev(rainbow(2000))[1001:2000]
  z_min <- min(triangles[,3])
  z_max <- max(triangles[,3])
  
  color_function <- function(v1, v2, v3){
    idx <- round( 1000 * ( (v3 - z_min)/(z_max - z_min)))
    return(color_vec[idx])
  }
  
  
  drawScene(makeTriangles(triangles, col.mesh = "black", color = color_function),
            R.mat = R_mat)
  
  
}
