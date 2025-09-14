###NULL MODEL 3 with GOWER DISTANCES AND PCOA SCORES----------------------------
##Simulates virtual communities by co-occurrence weighted drawing of species 
#from the native species pool.
#Beal's probability of species occurrence in the plot is used as a weight. 
#The higher the probability, the higher the chance of being selected.
#This function uses Gower's (1971) distance with Podani's (1999) extension for
#ordinal variables. Gower distance matrix is then ordinated using principal
#coordinate analysis and species scores on ordination axes are used instead of 
#original trait values.
#Returns distance statistics, their deviations from null expectation and 
#empirical p-values.

##Load packages-----------------------------------------------------------------
library(tidyverse)
library(ade4)
library(vegan)
library(reshape2)
library(parallel)
library(pbapply)
library(FD)
library(ape)

##Load functions----------------------------------------------------------------

source("./Functions/Edistance_weighted.R")
source("./Functions/stats.gpcoa.R")
source("./Functions/stats.gpcoa.calc.R")
source("./Functions/stats.gpcoa.ali.R")
source("./Functions/stats.gpcoa.ali.calc.R")
source("./Functions/null3.gpcoa.simulator.R")
source("./Functions/null3ali.gpcoa.simulator.R")

##Run simulation for native vs alien (naturalized and invasive) species--------- 

tr.vars <- c("HEIGHT", "SLA", "Seed.mass", "Leaf.area", "LDMC.mean", "FLOWERING.MEAN", "FLOWERING.LENGTH", "Gen.size.2C")

NULL3.gpcoa.missFor <- list(T=list(), X=list(), S=list(), M=list(), K=list(), L=list())

for(q in c("T", "X", "S", "M", "K", "L"))
{
  print(q)
  
  if(q %in% c("T", "X", "S", "M"))
  {
    print("Running simulations...")
    
    spe.nat <- CommData.missFor.ord %>% filter(Habitat == q & INVASION.STATUS == "native") %>%
      mutate(presence = 1) %>%
      reshape2::dcast(., ID ~ Species, value.var="presence", fill = 0) %>% 
      column_to_rownames(var="ID")
    
    b <- beals(spe.nat, include=F)*100
    bl <- SpecData.missFor.ord[colnames(b), "Life.form"]
    
    l <- split(CommData.missFor.ord[CommData.missFor.ord$Habitat == q,], CommData.missFor.ord[CommData.missFor.ord$Habitat == q, "ID"])
    invaded <- unlist(lapply(l, FUN = function(x){length(unique(x$NATIVE.ALIEN))}))
    l <- l[invaded == 2]
    
    #sqrt transform covers
    l <- lapply(l, FUN = function(x) {
      x$cover <- sqrt(x$cover)
      return(x)})
    
    cl <- makeCluster(10)
    clusterExport(cl, list("null3.gpcoa.simulator", "as.randtest", "stats.gpcoa", "stats.gpcoa.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames", "Edistance_weighted"))
    clusterSetRNGStream(cl,1234)
    
    d <- pblapply(l, FUN = null3.gpcoa.simulator, cl = cl, 
                  beals.scores = b, bl = bl, tr = SpecData.missFor.ord, tr.vars = tr.vars, nperm=999,
                  save.simulations = NULL)
    
    stopCluster(cl)
    
    NULL3.gpcoa.missFor[[q]] <- do.call(rbind.data.frame, d)
    
    rm(l)
  }
  
  if(q %in% c("K", "L"))
  {
    ht <- list(Herb=NULL, Tree=NULL)
    
    spe.nat <- CommData.missFor.ord %>% filter(Habitat == q & INVASION.STATUS == "native") %>%
      mutate(presence = 1) %>%
      reshape2::dcast(., ID ~ Species, value.var="presence", fill = 0) %>% 
      column_to_rownames(var="ID")
    
    b <- beals(spe.nat, include=F)*100
    bl <- SpecData.missFor.ord[colnames(b), "Life.form"]
    
    ##run simulations separately for herb and tree layer
    for(m in c("Herb", "Tree"))  
    {
      print(m)
      print("Running simulations...")
      
      l <- split(CommData.missFor.ord[CommData.missFor.ord$Habitat == q & CommData.missFor.ord$Herb.Tree == m,], CommData.missFor.ord[CommData.missFor.ord$Habitat == q & CommData.missFor.ord$Herb.Tree == m, "ID"])
      invaded <- unlist(lapply(l, FUN = function(x){length(unique(x$NATIVE.ALIEN))}))
      l <- l[invaded == 2]
      
      #sqrt transform covers
      l <- lapply(l, FUN = function(x) {
        x$cover <- sqrt(x$cover)
        return(x)})
      
      cl <- makeCluster(10)
      clusterExport(cl, list("null3.gpcoa.simulator", "as.randtest", "stats.gpcoa", "stats.gpcoa.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames", "Edistance_weighted"))
      clusterSetRNGStream(cl,1234)
      
      d <- pblapply(l, FUN = null3.gpcoa.simulator, cl = cl, 
                    beals.scores=b, bl=bl, tr=SpecData.missFor.ord, tr.vars = tr.vars, nperm=999,
                    save.simulations = NULL)
      
      stopCluster(cl)
      
      ht[[m]] <- do.call(rbind.data.frame, d)
      
      rm(l)
    }
    
    NULL3.gpcoa.missFor[[q]] <- ht
  }
  
}

saveRDS(NULL3.gpcoa.missFor, file="NULL3.gpcoa.missFor.rds")

###Run simulation for naturalized vs invasive species--------------------------- 

natz.id <- unique(CommData.missFor.ord$ID[CommData.missFor.ord$INVASION.STATUS == "naturalized"])
natz.plots <- CommData.missFor.ord[CommData.missFor.ord$ID %in% natz.id, ]

NULL3ali.gpcoa.missFor <- list(T=list(), X=list(), S=list(), M=list(), K=list(), L=list())

for(q in c("K"))#c("T", "X", "S", "M", "K", "L")
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
    bl <- SpecData.missFor.ord[colnames(b), "Life.form"]
    bs <- SpecData.missFor.ord[colnames(b), "INVASION.STATUS"]
    
    l <- split(natz.plots[natz.plots$Habitat == q,], natz.plots[natz.plots$Habitat == q, "ID"])
    
    #sqrt transform covers
    l <- lapply(l, FUN = function(x) {
      x$cover <- sqrt(x$cover)
      return(x)})
    
    cl <- makeCluster(10)
    clusterExport(cl, list("null3ali.gpcoa.simulator", "as.randtest", "stats.gpcoa.ali", "stats.gpcoa.ali.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames", "Edistance_weighted"))
    clusterSetRNGStream(cl,1234)
    
    d <- pblapply(l, FUN = null3ali.gpcoa.simulator, cl = cl, 
                  beals.scores = b, bl = bl, bs = bs, tr = SpecData.missFor.ord, tr.vars = tr.vars, nperm=999,
                  save.simulations = NULL)
    
    stopCluster(cl)
    
    NULL3ali.gpcoa.missFor[[q]] <- do.call(rbind.data.frame, d)
    NULL3ali.gpcoa.missFor[[q]] <- NULL3ali.gpcoa.missFor[[q]][ NULL3ali.gpcoa.missFor[[q]]$INVASION.STATUS != "native", ]
    
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
    bl <- SpecData.missFor.ord[colnames(b), "Life.form"]
    bs <- SpecData.missFor.ord[colnames(b), "INVASION.STATUS"]
    
    ##run simulations separately for herb and tree layer
    for(m in c("Herb", "Tree"))
    {
      print(m)
      print("Running simulations...")
      
      l <- split(natz.plots[natz.plots$Habitat == q & natz.plots$Herb.Tree == m,], natz.plots[natz.plots$Habitat == q & natz.plots$Herb.Tree == m, "ID"])
      
      #sqrt transform covers
      l <- lapply(l, FUN = function(x) {
        x$cover <- sqrt(x$cover)
        return(x)})
      
      cl <- makeCluster(10)
      clusterExport(cl, list("null3ali.gpcoa.simulator", "as.randtest", "stats.gpcoa.ali", "stats.gpcoa.ali.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames", "Edistance_weighted"))
      clusterSetRNGStream(cl,1234)
      
      d <- pblapply(l, FUN = null3ali.gpcoa.simulator, cl = cl, 
                    beals.scores = b, bl = bl, bs = bs, tr = SpecData.missFor.ord, tr.vars = tr.vars, nperm=999,
                    save.simulations = NULL)
      
      stopCluster(cl)
      
      ht[[m]] <- do.call(rbind.data.frame, d)
      ht[[m]] <- ht[[m]][ ht[[m]]$INVASION.STATUS != "native", ]
      
      rm(l)
    }
    
    NULL3ali.gpcoa.missFor[[q]] <- ht
  }
  
}

saveRDS(NULL3ali.gpcoa.missFor, file="NULL3ali.gpcoa.missFor.rds")