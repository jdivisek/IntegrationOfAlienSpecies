########################################################################
###       SPECIES DISTRIBUTION IN THE TRAIT SPACE OF EACH PLOT       ###
###                AAC quantification and simulation                 ### 
###                   by unconstrained null model                    ###
########################################################################

##1) Quantify the distribution of native and alien (either naturalized or invasive)
##species in the eight-dimensional functional trait space of each plot.
##2) Account for random expectation using unconstrained null model which randomly
##samples species from habitat species pool

##load packages
library(tidyverse)
library(ade4)
library(vegan)
library(reshape2)
library(parallel)
library(pbapply)
library(rlist)

###AAC calculator
accum <- function(dat, centr, st)
{
  d <- as.matrix(dist(rbind(centr, as.data.frame(dat)), method="euclidean"))[-1,1]
  if(max(d) != 0) {d <- d/max(d)}
  
  acc <- sapply(sort(unique(st)), FUN = function(s){
    1 - (sum(sapply(seq(0.01, 1, 0.01), FUN=function(x){sum(d[st == s] <= x)/length(d[st == s])}))*0.01)
  })
  
  return(acc)
}

##Assemble AAC values for native and alien species
accum.calc <- function(one.plot, tr.vars, nat.means=NULL)
{
  if(length(unique(one.plot$INVASION.STATUS)) > 1 & any(one.plot$INVASION.STATUS %in% c("native")))
  {
    if(is.null(nat.means)){nat.means <- apply(one.plot[one.plot$INVASION.STATUS == "native", tr.vars], 2, mean)}
    
    n <- as.data.frame(table(one.plot$INVASION.STATUS))
    colnames(n) <- c("INVASION.STATUS", "N")
    n$acc <- accum(dat = one.plot[,tr.vars], 
                   centr = nat.means,
                   st = one.plot$INVASION.STATUS)
    n <- n %>% mutate(ID = one.plot$ID[1], .before=1)
    
    return(n)
  }
}

##Null model - unconstrained drawing from native species pool
accum.simulator <- function(one.plot, tr, tr.vars, nperm=1000, save.simulations=NULL)
{
  if(length(unique(one.plot$INVASION.STATUS)) > 1 & any(one.plot$INVASION.STATUS %in% c("native")))
  {
    op <- split(one.plot, one.plot$Life.form2)##split based on life forms
    
    op.sim <- list()
    step <- 1
    repeat
    {
      op.sim[[step]] <- lapply(op, FUN = function(x)
      {
        w <- rownames(tr)[tr[,x$Habitat[1]] > 0 & tr$Life.form2 == x$Life.form2[1] & tr$INVASION.STATUS == "native"]
        
        size <- nrow(x)
        if(size <= length(w)) {sel <- sample(w, size = size, replace=F)}
        else {sel <- sample(w, size = size, replace=T)}
        
        x[, tr.vars] <- tr[sel, tr.vars]
        
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
    
    obs <- accum.calc(one.plot =  one.plot, tr.vars = tr.vars, nat.means = nat.means)
    obs[, tr.vars] <- apply(one.plot[, tr.vars], 2, FUN=function(x){ 
      accum(x, centr=mean(x[one.plot$INVASION.STATUS == "native"]), st = one.plot$INVASION.STATUS)} )
    
    ###list of simulations
    op.sim <- split(op.sim[,c("INVASION.STATUS", tr.vars)], op.sim$simulation)
    
    sim <- lapply(op.sim, FUN = function(x){
      s <- data.frame(acc = accum(dat=x[, tr.vars], centr=apply(x[x$INVASION.STATUS == "native", tr.vars], 2, mean), st=x$INVASION.STATUS))
      
      s[, tr.vars] <- apply(x[, tr.vars], 2, FUN=function(y){
        accum(y, centr=mean(y[x$INVASION.STATUS == "native"]), st =x$INVASION.STATUS)})
      
      s <- mutate(.data = s, INVASION.STATUS = rownames(s), .before = 1)
      
      return(s[order(rownames(s)), ])
    })
    sim <- do.call(rbind.data.frame, sim)
    
    merged <- rbind(obs[, c("INVASION.STATUS", "acc", tr.vars)], sim)
    
    for(st in obs$INVASION.STATUS)
    {
      obs[obs$INVASION.STATUS == st, paste0(c("acc", tr.vars), ".obs.minus.exp")] <- apply(merged[merged$INVASION.STATUS == st, -1], 2, FUN = function(x){ x[1] - mean(x[-1])}) 
      obs[obs$INVASION.STATUS == st, paste0(c("acc", tr.vars), ".ses")] <- apply(merged[merged$INVASION.STATUS == st, -1], 2, FUN = function(x){ (x[1] - mean(x[-1]))/sd(x[-1])}) 
      obs[obs$INVASION.STATUS == st, paste0(c("acc", tr.vars), ".p.value")] <- apply(merged[merged$INVASION.STATUS == st, -1], 2, FUN = function(x){ 
        as.randtest(sim = x[-1], obs = x[1], alter = "two-sided")$pvalue}) 
    }
    
    if(!is.null(save.simulations))
    {
      write.table(sim, file=paste0(save.simulations, "/", one.plot$ID[1], ".txt"), sep="\t", dec=".")
    }
    
    return(obs)
    
  }
}

###RUN SIMULATION--------------------------------------------------------------- 

ACC.NULL <- list(T=list(), X=list(), S=list(), M=list(), K=list(), L=list())

for(q in c("T", "X", "S", "M", "K", "L"))#
{
  print(q)
  
  if(q %in% c("T", "X", "S", "M"))
  {
    print("Running simulations...")
    
    l <- split(spel.imp.ls[spel.imp.ls$Habitat == q,], spel.imp.ls[spel.imp.ls$Habitat == q, "ID"])
    
    cl <- makeCluster(16)
    clusterExport(cl, list("accum.simulator", "as.randtest", "accum", "accum.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames"))
    clusterSetRNGStream(cl,1234)
    
    d <- pblapply(l, FUN = accum.simulator, cl = cl, 
                  tr=trait.imp.ls, tr.vars = tr.vars, nperm=1000,
                  save.simulations = NULL)
    
    stopCluster(cl)
    
    ACC.NULL[[q]] <- do.call(rbind.data.frame, d)
    
    rm(l)
  }
  
  if(q %in% c("K", "L"))
  {
    ht <- list(Herb=NULL, Tree=NULL)
    
    for(m in c("Herb", "Tree"))
    {
      print(m)
      print("Running simulations...")
      
      l <- split(spel.imp.ls[spel.imp.ls$Habitat == q & spel.imp.ls$Herb.Tree == m,], spel.imp.ls[spel.imp.ls$Habitat == q & spel.imp.ls$Herb.Tree == m, "ID"])
      
      cl <- makeCluster(16)
      clusterExport(cl, list("accum.simulator", "as.randtest", "accum", "accum.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames"))
      clusterSetRNGStream(cl,1234)
      
      d <- pblapply(l, FUN = accum.simulator, cl = cl, 
                    tr=trait.imp.ls, tr.vars = tr.vars, nperm=1000,
                    save.simulations = NULL)
      
      stopCluster(cl)
      
      ht[[m]] <- do.call(rbind.data.frame, d)
      
      rm(l)
    }
    
    ACC.NULL[[q]] <- ht
  }
}

