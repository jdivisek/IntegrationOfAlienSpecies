########################################################################
###       SPECIES DISTRIBUTION IN THE TRAIT SPACE OF EACH PLOT       ###
###          d and E-distance quantification and simulation          ### 
###                   by constrained null model                      ###
########################################################################

##1) Quantify the distribution of naturalized and invasive
##species in the eight-dimensional functional trait space of each plot.
##2) Account for "random" expectation using null model which samples species
##naturalized species from habitat species pool based on their co-occurrence
##with native species

##load packages
library(tidyverse)
library(ade4)
library(vegan)
library(reshape2)
library(parallel)
library(pbapply)
library(rlist)
library(energy)

###AAC calculator
stats.ali <- function(dat, centr, st)
{
  d <- as.matrix(dist(rbind(centr, as.data.frame(dat)), method="euclidean"))
  if(max(d) != 0) {d <- d/max(d)}
  
  stu <- sort(unique(st))
  
  ##Relativize distances from the center
  d1 <- d[-1,1]
  if(max(d1) != 0) {d1 <- d1/max(d1)}

  ##Mean distance from the center
  dm <- sapply(stu, FUN = function(s){mean(d1[st == s])})

  ##Make pairwise distances relative
  d <- d[-1,-1]
  if(max(d) != 0) {d <- d/max(d)}
  
  ##Energy distance
  d1 <- d[order(st), order(st)]
  Edist <- as.matrix(energy::edist(d1, table(st), distance=TRUE, alpha=1))
  dimnames(Edist) <- list(stu, stu)
  Edist <- Edist["naturalized", ]
  
  ##MPD - within groups
  mpd <- sapply(stu, FUN = function(s){
    if(sum(st == s) == 1){ return(NA) }
    else { mean(as.numeric(as.dist(d[st == s, st == s])))}})
  
  diag(d) <- NA
  ##MPD - between groups
  mpd2n <- sapply(stu, FUN = function(s){ 
    if(sum(st == "naturalized") == 1){ mean(as.matrix(d[st == "naturalized", st == s]), na.rm=T) }
    else { 
      mean(as.matrix(d[st == "naturalized", st == s]), na.rm=T) }})
  
  ##MNTD - within groups
  mntd <- sapply(stu, FUN = function(s){ 
    if(sum(st == s) == 1){ return(NA) }
    else { 
      mean(apply(as.data.frame(d[st == s, st == s]), 2, min, na.rm=T)) }})
  
  ##MNTD - between groups
  mnntd <- sapply(stu, FUN = function(s){ 
    if(sum(st == "naturalized") == 1){ mean(as.matrix(d[st == "naturalized", st == s]), na.rm=T) }
    else { 
      mean(apply(as.data.frame(d[st == "naturalized", st == s]), 2, min, na.rm=T)) }})
  
  return(cbind(dm, Edist, mnntd, mntd, mpd2n, mpd))
}


##Assemble distance values for native and alien species
stats.ali.calc <- function(one.plot, tr.vars, nat.means=NULL)
{
  if(length(unique(one.plot$INVASION.STATUS)) > 2 & any(one.plot$INVASION.STATUS %in% c("native", "naturalized")))
  {
    if(is.null(nat.means)){nat.means <- apply(one.plot[one.plot$INVASION.STATUS == "native", tr.vars], 2, mean)}
    
    n <- as.data.frame(table(one.plot$INVASION.STATUS))
    colnames(n) <- c("INVASION.STATUS", "N")
    n <- cbind(n, stats.ali(dat = one.plot[,tr.vars], 
                            centr = nat.means,
                            st = one.plot$INVASION.STATUS))
    n <- n %>% mutate(ID = one.plot$ID[1], .before=1)
    
    return(n)
  }
}

##Null model - constrained by species co-occurrence
##co-occurrence based drawing from naturalized species pool
stats.ali.w.simulator <- function(one.plot, beals.scores, bl, bs, tr, tr.vars, nperm=999, save.simulations=NULL)
{
  if(length(unique(one.plot$INVASION.STATUS)) > 2 & any(one.plot$INVASION.STATUS %in% c("native", "naturalized")))
  {
    op <- split(one.plot, one.plot$Life.form2)##split based on life forms
    
    op.sim <- list()
    step <- 1
    repeat
    {
      op.sim[[step]] <- lapply(op, FUN = function(x)
      {
        ###sampling for alien species
        if(any(x$NATIVE.ALIEN == "alien"))
        {
          w <- setNames(beals.scores[as.character(x$ID[1]), bl==x$Life.form2[1] & bs == "naturalized"],
                        colnames(beals.scores)[bl==x$Life.form2[1] & bs == "naturalized"])
          if(sum(w) > 0){
            w <- w[w > 0] #remove species with Beals score < 0
            # w <- w[setdiff(names(w), one.plot$Species[one.plot$INVASION.STATUS == "native" & one.plot$Life.form2 == x$Life.form2[1]] )]
          } #exclude species already in the plot
          else{ w[names(w)] <- 1 }
          
          size <- sum(x$INVASION.STATUS != "native")
          if(size <= length(w)) {sel <- sample(names(w), size = size, prob=w, replace=F)}
          else {sel <- sample(names(w), size = size, prob=w, replace=T)}
          
          x[x$INVASION.STATUS != "native", tr.vars] <- tr[sel, tr.vars]
          sel <- NULL
        }
        
        ###sampling for native species
        w <- setNames(beals.scores[as.character(x$ID[1]), bl==x$Life.form2[1] & bs == "native"],
                      colnames(beals.scores)[bl==x$Life.form2[1] & bs == "native"])
        if(sum(w) > 0){
          w <- w[w > 0] #remove species with Beals score < 0
          # w <- w[setdiff(names(w), one.plot$Species[one.plot$INVASION.STATUS == "native" & one.plot$Life.form2 == x$Life.form2[1]] )]
        } #exclude species already in the plot
        else{ w[names(w)] <- 1 }
        
        size <- sum(x$INVASION.STATUS == "native")
        if(size <= length(w)) {sel <- sample(names(w), size = size, prob=w, replace=F)}
        else {sel <- sample(names(w), size = size, prob=w, replace=T)}
        
        x[x$INVASION.STATUS == "native", tr.vars] <- tr[sel, tr.vars]
        
        return(x)})
      
      op.sim[[step]] <- do.call(rbind.data.frame, op.sim[[step]]) %>% mutate(simulation = step)
      
      step <- step+1
      if(step == nperm+1)
      {
        break
      }
      
    }
    
    op.sim <- do.call(rbind.data.frame, op.sim)
    
    nat.means <- apply(one.plot[one.plot$INVASION.STATUS == "native", tr.vars], 2, mean)
    
    obs <- stats.ali.calc(one.plot = one.plot, tr.vars = tr.vars, nat.means = nat.means)
    obs.t <- apply(one.plot[, tr.vars], 2, FUN=function(x){ 
      stats.ali(x, centr=mean(x[one.plot$INVASION.STATUS == "native"]), st = one.plot$INVASION.STATUS)}, simplify = F )
    obs <- cbind(obs, do.call(cbind.data.frame, obs.t))
    
    ###list of simulations
    op.sim <- split(op.sim[,c("INVASION.STATUS", tr.vars)], op.sim$simulation)
    
    sim <- lapply(op.sim, FUN = function(x){
      # stats.ali(dat=x[, tr.vars], centr=nat.means, st=x$INVASION.STATUS)[s]
      si <- stats.ali(dat=x[, tr.vars], centr=apply(x[x$INVASION.STATUS == "native", tr.vars], 2, mean), st=x$INVASION.STATUS)
      
      si.t <- apply(x[, tr.vars], 2, FUN=function(y){
        stats.ali(y, centr=mean(y[x$INVASION.STATUS == "native"]), st =x$INVASION.STATUS)}, simplify = F)
      si <- cbind(si, do.call(cbind.data.frame, si.t))
      
      si <- mutate(.data = si, INVASION.STATUS = rownames(si), .before = 1)
      
      return(si[order(rownames(si)), ])
    })

    sim <- do.call(rbind.data.frame, sim)
    
    merged <- rbind(obs[, !colnames(obs) %in% c("ID", "N")], sim)
    
    for(st in obs$INVASION.STATUS)
    {
      ## Observed minus expected values
      obs[obs$INVASION.STATUS == st, paste0(colnames(merged)[-1], "_O-E")] <- apply(merged[merged$INVASION.STATUS == st, -1], 2, FUN = function(x){ x[1] - mean(x[-1])}) 
      
      ##Standardized effect size
      # obs[obs$INVASION.STATUS == st, paste0(colnames(merged)[-1], "_SES")] <- apply(merged[merged$INVASION.STATUS == st, -1], 2, FUN = function(x){ (x[1] - mean(x[-1]))/sd(x[-1])}) 
      
      ##Empirical p-value for distances
      obs[obs$INVASION.STATUS == st, paste0(colnames(merged)[-1], "_P-val")] <- apply(merged[merged$INVASION.STATUS == st, -1], 2, FUN = function(x){ 
        as.randtest(sim = x[-1], obs = x[1], alter = "greater")$pvalue}) 
      
    }
    
    if(!is.null(save.simulations))
    {
      write.table(sim, file=paste0(save.simulations, "/", one.plot$ID[1], ".txt"), sep="\t", dec=".")
    }
    
    return(obs)
    
  }
  
}

###RUN SIMULATION--------------------------------------------------------------- 

STATali.NULLw <- list(T=list(), X=list(), S=list(), M=list(), K=list(), L=list())

natz.id <- unique(spel.imp.ls$ID[spel.imp.ls$INVASION.STATUS == "naturalized"])
natz.plots <- spel.imp.ls[spel.imp.ls$ID %in% natz.id, ]

for(q in c("T", "X", "S", "M","K", "L"))# 
{
  print(q)
  
  if(q %in% c("T", "X", "S", "M"))
  {
    print("Running simulations...")
    
    spe.nat <- natz.plots %>% filter(Habitat == q & INVASION.STATUS %in% c("native", "naturalized")) %>% 
      mutate(presence = 1) %>%
      reshape2::dcast(., ID ~ Species, value.var="presence", fill = 0) %>% 
      column_to_rownames(var="ID")
    
    b <- beals(spe.nat, include=F)*100
    bl <- trait.imp.ls[colnames(b), "Life.form2"]
    bs <- trait.imp.ls[colnames(b), "INVASION.STATUS"]
    
    l <- split(natz.plots[natz.plots$Habitat == q,], natz.plots[natz.plots$Habitat == q, "ID"])
    
    cl <- makeCluster(10)
    clusterExport(cl, list("stats.ali.w.simulator", "as.randtest", "stats.ali", "stats.ali.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames"))
    clusterSetRNGStream(cl,1234)
    
    d <- pblapply(l, FUN = stats.ali.w.simulator, cl = cl, 
                  beals.scores=b, bl=bl, bs=bs, tr=trait.imp.ls, tr.vars = tr.vars, nperm=999,
                  save.simulations = NULL)#
    
    stopCluster(cl)
    
    STATali.NULLw[[q]] <- do.call(rbind.data.frame, d)
    STATali.NULLw[[q]] <- STATali.NULLw[[q]][ STATali.NULLw[[q]]$INVASION.STATUS != "native", ]
    
    rm(l)
  }
  
  if(q %in% c("K", "L"))
  {
    ht <- list(Herb=NULL, Tree=NULL)
    
    spe.nat <- natz.plots %>% filter(Habitat == q & INVASION.STATUS %in% c("native", "naturalized")) %>%
      mutate(presence = 1) %>%
      reshape2::dcast(., ID ~ Species, value.var="presence", fill = 0) %>% 
      column_to_rownames(var="ID")
    
    b <- beals(spe.nat, include=F)*100
    bl <- trait.imp.ls[colnames(b), "Life.form2"]
    bs <- trait.imp.ls[colnames(b), "INVASION.STATUS"]
    
    for(m in c("Herb", "Tree"))
    {
      print(m)
      print("Running simulations...")
      
      l <- split(natz.plots[natz.plots$Habitat == q & natz.plots$Herb.Tree == m,], natz.plots[natz.plots$Habitat == q & natz.plots$Herb.Tree == m, "ID"])
      
      cl <- makeCluster(10)
      clusterExport(cl, list("stats.ali.w.simulator", "as.randtest", "stats.ali", "stats.ali.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames"))
      clusterSetRNGStream(cl,1234)
      
      d <- pblapply(l, FUN = stats.ali.w.simulator, cl = cl, 
                    beals.scores=b, bl=bl, bs=bs, tr=trait.imp.ls, tr.vars = tr.vars, nperm=999,
                    save.simulations = NULL)
      
      stopCluster(cl)
      
      ht[[m]] <- do.call(rbind.data.frame, d)
      ht[[m]] <- ht[[m]][ ht[[m]]$INVASION.STATUS != "native", ]
      
      rm(l)
    }
    
    STATali.NULLw[[q]] <- ht
  }
  
  save.image("~/R Working Directory/Trait simulations 2024/.RData")
}
