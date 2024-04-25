########################################################################
###       SPECIES DISTRIBUTION IN THE TRAIT SPACE OF EACH PLOT       ###
###        AAC quantification for equalized number of species        ###
###             and simulation by constrained null model             ###
########################################################################

##1) Quantify the distribution of native and alien (either naturalized or invasive)
##species in the eight-dimensional functional trait space of each plot - EQUAL
##NUMBER OF NATIVE AND ALIEN SPECIES IN CONSIDERED IN EACH PLOT - ONLY MULTIDIMENSIONAL 
##AAC VALUES ARE CALCULATED
##2) Account for "random" expectation using null model which samples species from
##habitat species pool based on their co-occurrence

##load packages
library(tidyverse)
library(ade4)
library(vegan)
library(reshape2)
library(parallel)
library(pbapply)
library(rlist)

##This code samples equal number of native and alien species in a community
##It can handle ONLY TWO INVASION STATUSES
sample.eq <- function(y)
{
  size <- min(table(y$INVASION.STATUS))
  y %>% group_by(INVASION.STATUS) %>% sample_n(size = size) %>% ungroup() %>% as.data.frame() %>% return()
}

###AAC calculator
accum.eq <- function(dat, centr, d.max, st)
{
  d <- as.matrix(dist(rbind(centr, as.data.frame(dat)), method="euclidean"))[-1,1]
  if(max(d) != 0) {d <- d/d.max}
  
  acc <- sapply(sort(unique(st)), FUN = function(s){
    1 - (sum(sapply(seq(0.01, 1, 0.01), FUN=function(x){sum(d[st == s] <= x)/length(d[st == s])}))*0.01)
  })
  
  return(acc)
}

##Assemble AAC values for native and alien species
accum.eq.calc <- function(one.plot, tr.vars, nat.means=NULL, d.max)
{
  if(length(unique(one.plot$INVASION.STATUS)) > 1 & any(one.plot$INVASION.STATUS %in% c("native")))
  {
    n <- as.data.frame(table(one.plot$INVASION.STATUS))
    colnames(n) <- c("INVASION.STATUS", "N")
    n$acc <- accum.eq(dat = one.plot[,tr.vars], 
                      centr = nat.means,
                      st = one.plot$INVASION.STATUS,
                      d.max = d.max)
    n <- n %>% mutate(ID = one.plot$ID[1], .before=1)
    
    return(n)
  }
}

##Null model - constrained by species co-occurrence
##co-occurrence based drawing from native species pool
accum.w.eq.simulator <- function(one.plot, beals.scores, bl, tr, tr.vars, nperm=1000, save.simulations=NULL)
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
        w <- setNames(beals.scores[as.character(x$ID[1]), bl==x$Life.form2[1]],
                      colnames(beals.scores)[bl==x$Life.form2[1]])
        if(sum(w) > 0){
          w <- w[w > 0] #remove species with Beals score < 0
          # w <- w[setdiff(names(w), one.plot$Species[one.plot$INVASION.STATUS == "native" & one.plot$Life.form2 == x$Life.form2[1]] )]
        } #exclude species already in the plot
        else{ w[names(w)] <- 1 }
        
        size <- nrow(x)
        if(size <= length(w)) {sel <- sample(names(w), size = size, prob=w, replace=F)}
        else {sel <- sample(names(w), size = size, prob=w, replace=T)}
        
        x[, tr.vars] <- tr[sel, tr.vars]
        
        return(x)})
      
      op.sim[[step]] <- do.call(rbind.data.frame, op.sim[[step]]) %>% mutate(simulation = step)
      
      step <- step+1
      if(step == nperm+1)
      {
        break
      }
      
    }
    rm(op)
    
    op.sim <- do.call(rbind.data.frame, op.sim)
    
    nat.means <- apply(one.plot[one.plot$INVASION.STATUS == "native", tr.vars], 2, mean)
    d.max <- max(dist(rbind(nat.means, one.plot[, tr.vars]), method="euclidean"))
    
    ###list of simulations
    op.sim <- split(op.sim[, c("Species", "INVASION.STATUS", tr.vars)], op.sim$simulation)
    op.sim <- lapply(op.sim, FUN=function(x){
      x %>% remove_rownames() %>% column_to_rownames(var = "Species") %>% return()
    })
    
    if(any(one.plot$INVASION.STATUS == "invasive"))
    {
      sim.inv <- lapply(op.sim, FUN = function(x){
        op <- sample.eq(one.plot[one.plot$INVASION.STATUS %in% c("native", "invasive"), ])
        obs <- accum.eq.calc(one.plot = op, tr.vars = tr.vars, nat.means = nat.means, d.max=d.max)
        
        nat.means.sim <- apply(x[x$INVASION.STATUS == "native", tr.vars], 2, mean)
        d.max.sim <- max(dist(rbind(nat.means.sim, x[, tr.vars]), method="euclidean"))  
        
        obs$sim <- accum.eq(dat=x[op$Species, tr.vars], centr=nat.means.sim, d.max=d.max.sim, st=x[op$Species, "INVASION.STATUS"])
        obs <- as.data.frame(table(one.plot$INVASION.STATUS)) %>% left_join(x=obs, y=., by = join_by(INVASION.STATUS == Var1))
        obs$key <- "inv"
        
        return(obs)
      })
    }
    else{sim.inv <- NULL}
    
    if(any(one.plot$INVASION.STATUS == "naturalized"))
    {
      sim.natz <- lapply(op.sim, FUN = function(x){
        op <- sample.eq(one.plot[one.plot$INVASION.STATUS %in% c("native", "naturalized"), ])
        obs <- accum.eq.calc(one.plot = op, tr.vars = tr.vars, nat.means = nat.means, d.max=d.max)
        
        nat.means.sim <- apply(x[x$INVASION.STATUS == "native", tr.vars], 2, mean)
        d.max.sim <- max(dist(rbind(nat.means.sim, x[, tr.vars]), method="euclidean"))  
        
        obs$sim <- accum.eq(dat=x[op$Species, tr.vars], centr=nat.means.sim, d.max=d.max.sim, st=x[op$Species, "INVASION.STATUS"])
        obs <- as.data.frame(table(one.plot$INVASION.STATUS)) %>% left_join(x=obs, y=., by = join_by(INVASION.STATUS == Var1))
        obs$key <- "natz"
        
        return(obs)
      })
    }
    else{sim.natz <- NULL}
    
    sim <- append(sim.inv, sim.natz)
    sim <- do.call(rbind.data.frame, sim)
    
    sim$Diff <- sim$acc - sim$sim ##observed minus expected
    
    std.dev <- aggregate(sim ~ ID + INVASION.STATUS + key, data = sim, sd)
    
    res <- aggregate(cbind(N, acc, sim, Freq, Diff) ~ ID + INVASION.STATUS + key, data = sim, mean)
    
    res$SES <- res$Diff / std.dev$sim
    
    
    if(!is.null(save.simulations))
    {
      write.table(sim, file=paste0(save.simulations, "/", one.plot$ID[1], ".txt"), sep="\t", dec=".")
    }
    
    return(res)
    
  }
  
}

###RUN SIMULATION--------------------------------------------------------------- 

ACC.NULLw.eq <- list(T=list(), X=list(), S=list(), M=list(), K=list(), L=list())

for(q in c("T", "X", "S", "M", "K", "L"))#
{
  print(q)
  
  if(q %in% c("T", "X", "S", "M"))
  {
    print("Running simulations...")
    
    spe.nat <- spel.imp.ls %>% filter(Habitat == q & INVASION.STATUS == "native") %>%
      mutate(presence = 1) %>%
      reshape2::dcast(., ID ~ Species, value.var="presence", fill = 0) %>% 
      column_to_rownames(var="ID")
    
    b <- beals(spe.nat, include=F)*100
    bl <- trait.imp.ls[colnames(b), "Life.form2"]
    
    l <- split(spel.imp.ls[spel.imp.ls$Habitat == q,], spel.imp.ls[spel.imp.ls$Habitat == q, "ID"])
    
    cl <- makeCluster(16)
    clusterExport(cl, list("accum.w.eq.simulator", "accum.eq", "accum.eq.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames", "group_by", "sample_n", "sample.eq", "ungroup", "left_join", "join_by"))
    clusterSetRNGStream(cl,1234)
    
    d <- pblapply(l, FUN = accum.w.eq.simulator, cl = cl, 
                  beals.scores=b, bl=bl, tr=trait.imp.ls, tr.vars = tr.vars, nperm=1000,
                  save.simulations = NULL)
    
    stopCluster(cl)
    
    ACC.NULLw.eq[[q]] <- do.call(rbind.data.frame, d)
    
    rm(l)
  }
  
  if(q %in% c("K", "L"))
  {
    ht <- list(Herb=NULL, Tree=NULL)
    
    spe.nat <- spel.imp.ls %>% filter(Habitat == q & INVASION.STATUS == "native") %>%
      mutate(presence = 1) %>%
      reshape2::dcast(., ID ~ Species, value.var="presence", fill = 0) %>% 
      column_to_rownames(var="ID")
    
    b <- beals(spe.nat, include=F)*100
    bl <- trait.imp.ls[colnames(b), "Life.form2"]
    
    for(m in c("Herb", "Tree"))
    {
      print(m)
      print("Running simulations...")
      
      l <- split(spel.imp.ls[spel.imp.ls$Habitat == q & spel.imp.ls$Herb.Tree == m,], spel.imp.ls[spel.imp.ls$Habitat == q & spel.imp.ls$Herb.Tree == m, "ID"])
      
      cl <- makeCluster(16)
      clusterExport(cl, list("accum.w.eq.simulator", "accum.eq", "accum.eq.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames", "group_by", "sample_n", "sample.eq", "ungroup", "left_join", "join_by"))
      clusterSetRNGStream(cl,1234)
      
      d <- pblapply(l, FUN = accum.w.eq.simulator, cl = cl, 
                    beals.scores=b, bl=bl, tr=trait.imp.ls, tr.vars = tr.vars, nperm=1000,
                    save.simulations = NULL)
      
      stopCluster(cl)
      
      ht[[m]] <- do.call(rbind.data.frame, d)
      
      rm(l)
    }
    
    ACC.NULLw.eq[[q]] <- ht
  }
}

