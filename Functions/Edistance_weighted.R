##Returns unweighted or weighted E-distances (energy statistics) between species groups.

#@test.group = First of the compared groups (native, naturalized or invasive).
#@ref.group = Second of the compared groups (native, naturalized or invasive).
#@D = Pairwise distance matrix.
#@st = Character vector of status labels (i.e. native, naturalized and invasive).
#@covers = Numeric vector of species covers/abundances.
#@weighted = If TRUE, the cover-weighted version of the E-distance is returned.
#@coef = For comparability with energy::edist function. Coefficient for weighting the E-distance. See R documentation for the edist function. If "none" (default), no coefficient is used. 

Edistance_weighted <- function(
    test.group, ref.group, D, st, covers = NULL, weighted = FALSE,
    coef = c("none", "cluster", "disco", "mean.cover", "sum.cover")) 
{
  if (inherits(D, "dist")) {
    D <- as.matrix(D)
  }
  
  stopifnot(length(st) == nrow(D))
  coef <- match.arg(coef)
  
  N <- length(st)
  n1 <- sum(st == ref.group)
  n2 <- sum(st == test.group)
  st <- as.character(st)
  X_idx <- which(st == ref.group)
  Y_idx <- which(st == test.group)
  
  A <- mean(D[X_idx, Y_idx])
  B <- mean(D[X_idx, X_idx])
  C <- mean(D[Y_idx, Y_idx])
  
  if(is.null(covers) | weighted == FALSE)
  {
    E <- 2*A - B - C
    
  } else {
    
    W <- outer(covers, covers, FUN = "*")
    
    stopifnot(nrow(W) == nrow(D))
    
    WA <- weighted.mean(D[X_idx, Y_idx], w = W[X_idx, Y_idx])
    WB <- weighted.mean(D[X_idx, X_idx], w = W[X_idx, X_idx])
    WC <- weighted.mean(D[Y_idx, Y_idx], w = W[Y_idx, Y_idx])
    
    E <- 2*WA - WB - WC
  }
  
  if(!is.null(covers) & coef %in% c("mean.cover", "sum.cover"))
  {
    stopifnot(length(covers) == nrow(D))
    
    wx <- covers[X_idx]
    wy <- covers[Y_idx]
  }
  
  coefficient <- switch(
    coef,
    "none" = 1,
    "cluster" = (n1*n2)/(n1+n2),
    "disco" = (n1*n2)/(2*N),
    "mean.cover" = ifelse(!is.null(covers), ((mean(wx)*mean(wy))/(mean(wx)+mean(wy))), NA),
    "sum.cover" = ifelse(!is.null(covers), ((sum(wx)*sum(wy))/(sum(wx)+sum(wy))), NA)
  )
  
  result <- coefficient * E
  
  return(result)
}
