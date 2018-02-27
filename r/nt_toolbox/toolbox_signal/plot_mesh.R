#options(repr.plot.width=5, repr.plot.height=5)
options(warn=-1) # turns off warnings, to turn on: "options(warn=0)"

plot_mesh <- function (X, F, col="black"){
  ###
  
  ### plots the mesh defined by X : vertices and F : faces
  
  ###
  trimesh(t(F+1),t(X), col=col) 

}