
compute_huffcode <- function (pvecin) { 
  # Remove all zeros from vec - no code is needed for them 
  pvec=pvecin[which(pvecin>0)]; 
  
  # Make couple of checks 
  if (class(pvec) != "numeric") stop("Input must be numeric"); 
  if (length(pvec) < 2) stop("Numbers of positive elements must be 2 or more"); 
  
  # Recursion stop condition if the number of rows exactly 2 
  if (length(pvec) == 2) return(c("0","1")); 
  
  # Compress low probabilities values and call again. Use last min 
  min1loc <- length(pvec)-which.min(rev(pvec))+1; 
  min1 <- pvec[min1loc]; 
  pvec <- pvec[-min1loc]; 
  # min2loc <- which.min(pvec); 
  min2loc <- length(pvec)-which.min(rev(pvec))+1; 
  min2 <- pvec[min2loc]; 
  pvec <- pvec[-min2loc]; 
  pvec <- append(pvec,min1+min2); 
  
  code <- compute_huffcode(pvec); 
  
  # Get the last entry code, remove it and add split it back to previous entries 
  min1code <- paste0(tail(code,n=1),"0"); 
  min2code <- paste0(tail(code,n=1),"1"); 
  # Remove last entry and add this code to the min1, min2 locations 
  code <- code[-length(code)]; 
  # Add back the new codes 
  code <- append(code,min2code,min2loc-1); 
  code <- append(code,min1code,min1loc-1); 
  
  return(code); 
  
} 
