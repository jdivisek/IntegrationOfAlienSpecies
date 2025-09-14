##Returns distance statistics for native and alien (either naturalized or invasive) species in 
#the selected community (vegetation plot).

#@dat = Trait data (columns) for the species present in each community (rows).
#@st = Character vector of status labels (i.e. native, naturalized and invasive).
#@cover = Numeric vector of species covers/abundances.

stats <- function(dat, st, cover)
{
  dat <- as.data.frame(dat)
  
  ##calculate center of the trait space
  centr <- apply(dat[st == "native", , drop = FALSE], 2, FUN = mean)
  centr.w <- apply(dat[st == "native", , drop = FALSE], 2, FUN = weighted.mean, w = cover[st == "native"])
  
  d <- as.matrix(dist(rbind(centr, as.data.frame(dat)), method="euclidean")) #distances from unweighted mean (center of the trait space)
  dw <- as.matrix(dist(rbind(centr.w, as.data.frame(dat)), method="euclidean")) #distances from weighted mean (center of the trait space)
  
  stu <- sort(unique(st))
  
  ##make distances relative
  d1 <- d[-1,1]
  dw1 <- dw[-1,1]
  
  if(max(d1) != 0) {d1 <- d1/max(d1)}
  if(max(dw1) != 0) {dw1 <- dw1/max(dw1)}
  
  ##mean distance (D)
  dm <- sapply(stu, FUN = function(s){mean(d1[st == s])}) ##unweighted mean distance from unweighted centroid
  dwm <- sapply(stu, FUN = function(s){mean(dw1[st == s])}) ##unweighted mean distance from weighted centroid
  wdm <- sapply(stu, FUN = function(s){weighted.mean(d1[st == s], w = cover[st == s])}) #weighted mean distance from unweighted centroid
  wdwm <- sapply(stu, FUN = function(s){weighted.mean(dw1[st == s], w = cover[st == s])}) #weighted mean distance from weighted centroid
  
  #make distances relative
  d <- d[-1,-1]
  if(max(d) != 0) {d <- d/max(d)}
  
  ##Energy distance (E)
  Edist <- sapply(stu, FUN = Edistance_weighted, ref.group = "native", D = d, st = st, 
                  covers = NULL, weighted = FALSE, coef = "none")
  
  Edist.cluster <- sapply(stu, FUN = Edistance_weighted, ref.group = "native", D = d, st = st, 
                          covers = NULL, weighted = FALSE, coef = "cluster")
  
  Edist.mean.cover <- sapply(stu, FUN = Edistance_weighted, ref.group = "native", D = d, st = st, 
                             covers = cover, weighted = FALSE, coef = "mean.cover")
  
  Edist.sum.cover <- sapply(stu, FUN = Edistance_weighted, ref.group = "native", D = d, st = st, 
                            covers = cover, weighted = FALSE, coef = "sum.cover")
  
  wEdist <- sapply(stu, FUN = Edistance_weighted, ref.group = "native", D = d, st = st, 
                   covers = cover, weighted = TRUE, coef = "none")
  
  wEdist.cluster <- sapply(stu, FUN = Edistance_weighted, ref.group = "native", D = d, st = st, 
                           covers = cover, weighted = TRUE, coef = "cluster")
  
  wEdist.mean.cover <- sapply(stu, FUN = Edistance_weighted, ref.group = "native", D = d, st = st, 
                              covers = cover, weighted = TRUE, coef = "mean.cover")
  
  wEdist.sum.cover <- sapply(stu, FUN = Edistance_weighted, ref.group = "native", D = d, st = st, 
                             covers = cover, weighted = TRUE, coef = "sum.cover")
  
  return(cbind(dm, dwm, wdm, wdwm, Edist, Edist.cluster, Edist.mean.cover, Edist.sum.cover, wEdist, wEdist.cluster, wEdist.mean.cover, wEdist.sum.cover))
}
