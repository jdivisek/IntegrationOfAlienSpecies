##Returns distance statistics for naturalized and invasive species in 
#the selected community (vegetation plot).
#This function uses Gower's (1971) distance with Podani's (1999) extension for
#ordinal variables. Gower distance matrix is then ordinated using principal
#coordinate analysis and species scores on ordination axes are used instead of 
#original trait values.

#@dat = Trait data (columns) for the species present in each community (rows).
#@st = Character vector of status labels (i.e. native, naturalized and invasive).
#@cover = Numeric vector of species covers/abundances.

stats.gpcoa.ali <- function(dat, st, cover)
{
  require(FD)
  require(ape)
  
  dat <- as.data.frame(dat)
  
  ##check equal traits
  v <- apply(X = dat, MARGIN = 2, FUN = function(y){length(unique(y))})
  
  if(all(v > 1)) { 
    g <- FD::gowdis(x = dat, ord = "podani")
    
    ##pcoa of mixed-type traits
    p <- ape::pcoa(D = g, correction = "lingoes")
    
    if("vectors.cor" %in% names(p)){p <- p$vectors.cor
    } else {p <- p$vectors}
    
    centr <- apply(X = p[st == "native", , drop = FALSE], MARGIN = 2, FUN = mean)
    centr.w <- apply(X = p[st == "native", , drop = FALSE], MARGIN = 2, FUN = weighted.mean, w = cover[st == "native"])
    
    d <- as.matrix(dist(rbind(centr, as.data.frame(p))))
    dw <- as.matrix(dist(rbind(centr.w, as.data.frame(p))))
    
  } else {
    if(ncol(dat) > 1 & sum(v == 1) < ncol(dat)){
      dat <- dat[,v > 1]
      
      g <- FD::gowdis(x = dat, ord = "podani")
      
      ##pcoa of mixed-type traits
      p <- ape::pcoa(D = g, correction = "lingoes")
      
      if("vectors.cor" %in% names(p)){p <- p$vectors.cor
      } else {p <- p$vectors}
      
      centr <- apply(X = p[st == "native", , drop = FALSE], MARGIN = 2, FUN = mean)
      centr.w <- apply(X = p[st == "native", , drop = FALSE], MARGIN = 2, FUN = weighted.mean, w = cover[st == "native"])
      
      d <- as.matrix(dist(rbind(centr, as.data.frame(p))))
      dw <- as.matrix(dist(rbind(centr.w, as.data.frame(p))))
      
    } else {
      dat <- apply(X = dat, MARGIN = 2, FUN = as.numeric)
      
      centr <- apply(X = dat[st == "native", , drop = FALSE], MARGIN = 2, FUN = mean)
      centr.w <- apply(X = dat[st == "native", , drop = FALSE], MARGIN = 2, FUN = weighted.mean, w = cover[st == "native"])
      
      d <- as.matrix(dist(rbind(centr, as.data.frame(dat))))
      dw <- as.matrix(dist(rbind(centr.w, as.data.frame(dat))))
    }
  }
  
  stu <- sort(unique(st))
  
  ##make distances relative
  d1 <- d[-1,1]
  dw1 <- dw[-1,1]
  
  if(max(d1) != 0) {d1 <- d1/max(d1)}
  if(max(dw1) != 0) {dw1 <- dw1/max(dw1)}
  
  ##Mean distance (D)
  dm <- sapply(stu, FUN = function(s){mean(d1[st == s])}) ##unweighted mean distance from unweighted centroid
  dwm <- sapply(stu, FUN = function(s){mean(dw1[st == s])}) ##unweighted mean distance from weighted centroid
  wdm <- sapply(stu, FUN = function(s){weighted.mean(d1[st == s], w = cover[st == s])}) #weighted mean distance from unweighted centroid
  wdwm <- sapply(stu, FUN = function(s){weighted.mean(dw1[st == s], w = cover[st == s])}) #weighted mean distance from weighted centroid
  
  #make distances relative
  d <- d[-1,-1]
  if(max(d) != 0) {d <- d/max(d)}
  
  ##Energy distance (E)
  Edist <- sapply(stu, FUN = Edistance_weighted, ref.group = "naturalized", D = d, st = st, 
                  covers = NULL, weighted = FALSE, coef = "none")
  
  Edist.cluster <- sapply(stu, FUN = Edistance_weighted, ref.group = "naturalized", D = d, st = st, 
                          covers = NULL, weighted = FALSE, coef = "cluster")
  
  Edist.mean.cover <- sapply(stu, FUN = Edistance_weighted, ref.group = "naturalized", D = d, st = st, 
                             covers = cover, weighted = FALSE, coef = "mean.cover")
  
  Edist.sum.cover <- sapply(stu, FUN = Edistance_weighted, ref.group = "naturalized", D = d, st = st, 
                            covers = cover, weighted = FALSE, coef = "sum.cover")
  
  wEdist <- sapply(stu, FUN = Edistance_weighted, ref.group = "naturalized", D = d, st = st, 
                   covers = cover, weighted = TRUE, coef = "none")
  
  wEdist.cluster <- sapply(stu, FUN = Edistance_weighted, ref.group = "naturalized", D = d, st = st, 
                           covers = cover, weighted = TRUE, coef = "cluster")
  
  wEdist.mean.cover <- sapply(stu, FUN = Edistance_weighted, ref.group = "naturalized", D = d, st = st, 
                              covers = cover, weighted = TRUE, coef = "mean.cover")
  
  wEdist.sum.cover <- sapply(stu, FUN = Edistance_weighted, ref.group = "naturalized", D = d, st = st, 
                             covers = cover, weighted = TRUE, coef = "sum.cover")
  
  return(cbind(dm, dwm, wdm, wdwm, Edist, Edist.cluster, Edist.mean.cover, Edist.sum.cover, wEdist, wEdist.cluster, wEdist.mean.cover, wEdist.sum.cover))
}
