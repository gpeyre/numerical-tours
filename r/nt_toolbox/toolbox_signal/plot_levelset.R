



plot_levelset <- function(Z, f=c(), title="", lw=1.5, sbpt=c()){
  ####
  # f is supposed to be of the same shape as Z
  ####
  
  if (length(sbpt) >0){
    if (sbpt[3]==1){
      par(mfrow=sbpt[1:2]) 
    }
  }
  
  if (length(f)==0){
    f <- Z }

  Z <- as.cimg(t(Z))
  ct <- contours(Z, nlevels=1)
  imageplot(f, title, sbpt)
  purrr::walk(ct,function(v) lines(v$x,v$y,col="red",lw=lw))
}