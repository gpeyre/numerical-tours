library("colorRamps")
options(repr.plot.width=8, repr.plot.height=2.5)
options(warn=-1) # turns off warnings, to turn on: "options(warn=0)"

plot_spectogram <- function(S,title){
  S <- abs(S[1:(dim(S)[1]/2),])
  S <- log(S + 1e-4)
  image(t(S),col=matlab.like2(50), main=title)
}
