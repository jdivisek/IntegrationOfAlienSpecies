###NULL MODEL 3-----------------------------------------------------------------
##Simulates virtual communities by co-occurrence weighted drawing of species from 
#the native/naturalized species pool.
#Beal's probability of species occurrence in the plot is used as a weight. 
#The higher the probability, the higher the chance of being selected.
#Returns distance statistics, their deviations from null expectation and 
#empirical p-values.

##load packages-----------------------------------------------------------------
library(tidyverse)
library(ade4)
library(vegan)
library(reshape2)
library(parallel)
library(pbapply)

##Load functions----------------------------------------------------------------

source("./Functions/Edistance_weighted.R")
source("./Functions/stats.R")
source("./Functions/stats.calc.R")
source("./Functions/stats.ali.R")
source("./Functions/stats.ali.calc.R")
source("./Functions/null3.simulator.R")
source("./Functions/null3ali.simulator.R")

##Run simulation for native vs alien (naturalized and invasive) species--------- 

tr.vars <- c("HEIGHT", "SLA", "Seed.mass", "Leaf.area", "LDMC.mean", "FLOWERING.MEAN", "FLOWERING.LENGTH", "Gen.size.2C")

NULL3.missFor <- list(T=list(), X=list(), S=list(), M=list(), K=list(), L=list())

for(q in c("T", "X", "S", "M","K", "L"))
{
  print(q)

  if(q %in% c("T", "X", "S", "M"))
  {
    print("Running simulations...")

    #plot x native species matrix 
    spe.nat <- CommData.missFor %>% filter(Habitat == q & INVASION.STATUS == "native") %>%
      mutate(presence = 1) %>%
      reshape2::dcast(., ID ~ Species, value.var="presence", fill = 0) %>%
      column_to_rownames(var="ID")

    #beals co-occurrence scores
    b <- vegan::beals(spe.nat, include=F)*100
    bl <- SpecData.missFor[colnames(b), "Life.form"]

    l <- split(CommData.missFor[CommData.missFor$Habitat == q,], CommData.missFor[CommData.missFor$Habitat == q, "ID"])
    invaded <- unlist(lapply(l, FUN = function(x){length(unique(x$NATIVE.ALIEN))}))
    l <- l[invaded == 2]

    #sqrt transform covers
    l <- lapply(l, FUN = function(x) {
      x$cover <- sqrt(x$cover)
      return(x)})

    cl <- makeCluster(10)
    clusterExport(cl, list("null3.simulator", "as.randtest", "stats", "stats.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames", "Edistance_weighted"))
    clusterSetRNGStream(cl,1234)

    d <- pblapply(l, FUN = null3.simulator, cl = cl,
                  beals.scores = b, bl = bl, tr = SpecData.missFor, tr.vars = tr.vars, nperm=999,
                  save.simulations = NULL)

    stopCluster(cl)

    NULL3.missFor[[q]] <- do.call(rbind.data.frame, d)

    rm(l)
  }

  if(q %in% c("K", "L"))
  {
    ht <- list(Herb=NULL, Tree=NULL)

    #plot x native species matrix 
    spe.nat <- CommData.missFor %>% filter(Habitat == q & INVASION.STATUS == "native") %>%
      mutate(presence = 1) %>%
      reshape2::dcast(., ID ~ Species, value.var="presence", fill = 0) %>%
      column_to_rownames(var="ID")

    #beals co-occurrence scores
    b <- vegan::beals(spe.nat, include=F)*100
    bl <- SpecData.missFor[colnames(b), "Life.form"]

    ##run simulations separately for herb and tree layer
    for(m in c("Herb", "Tree"))
    {
      print(m)
      print("Running simulations...")

      l <- split(CommData.missFor[CommData.missFor$Habitat == q & CommData.missFor$Herb.Tree == m,], CommData.missFor[CommData.missFor$Habitat == q & CommData.missFor$Herb.Tree == m, "ID"])
      invaded <- unlist(lapply(l, FUN = function(x){length(unique(x$NATIVE.ALIEN))}))
      l <- l[invaded == 2]

      #sqrt transform covers
      l <- lapply(l, FUN = function(x) {
        x$cover <- sqrt(x$cover)
        return(x)})

      cl <- makeCluster(10)
      clusterExport(cl, list("null3.simulator", "apply", "as.randtest", "stats", "stats.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames", "Edistance_weighted"))
      clusterSetRNGStream(cl,1234)

      d <- pblapply(l, FUN = null3.simulator, cl = cl,
                    beals.scores = b, bl = bl, tr = SpecData.missFor, tr.vars = tr.vars, nperm=999,
                    save.simulations = NULL)

      stopCluster(cl)

      ht[[m]] <- do.call(rbind.data.frame, d)

      rm(l)
    }

    NULL3.missFor[[q]] <- ht
  }

}

saveRDS(NULL3.missFor, file="NULL3.missFor.rds")

###Run simulation for naturalized vs invasive species--------------------------- 

NULL3ali.missFor <- list(T=list(), X=list(), S=list(), M=list(), K=list(), L=list())

natz.id <- unique(CommData.missFor$ID[CommData.missFor$INVASION.STATUS == "naturalized"])
natz.plots <- CommData.missFor[CommData.missFor$ID %in% natz.id, ]

for(q in c("T", "X", "S", "M","K", "L"))
{
  print(q)
  
  if(q %in% c("T", "X", "S", "M"))
  {
    print("Running simulations...")
    
    #plot x native & naturalized species matrix 
    spe.nat <- natz.plots %>% filter(Habitat == q & INVASION.STATUS %in% c("native", "naturalized")) %>% 
      mutate(presence = 1) %>%
      reshape2::dcast(., ID ~ Species, value.var="presence", fill = 0) %>% 
      column_to_rownames(var="ID")
    
    #beals co-occurrence scores
    b <- vegan::beals(spe.nat, include=F)*100
    bl <- SpecData.missFor[colnames(b), "Life.form"]
    bs <- SpecData.missFor[colnames(b), "INVASION.STATUS"]
    
    l <- split(natz.plots[natz.plots$Habitat == q,], natz.plots[natz.plots$Habitat == q, "ID"])
    
    #sqrt transform covers
    l <- lapply(l, FUN = function(x) {
      x$cover <- sqrt(x$cover)
      return(x)})
    
    cl <- makeCluster(10)
    clusterExport(cl, list("null3ali.simulator", "as.randtest", "stats.ali", "stats.ali.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames", "Edistance_weighted"))
    clusterSetRNGStream(cl,1234)
    
    d <- pblapply(l, FUN = null3ali.simulator, cl = cl, 
                  beals.scores = b, bl = bl, bs = bs, tr = SpecData.missFor, tr.vars = tr.vars, nperm=999,
                  save.simulations = NULL)
    
    stopCluster(cl)
    
    NULL3ali.missFor[[q]] <- do.call(rbind.data.frame, d)
    NULL3ali.missFor[[q]] <- NULL3ali.missFor[[q]][ NULL3ali.missFor[[q]]$INVASION.STATUS != "native", ]
    
    rm(l)
  }
  
  if(q %in% c("K", "L"))
  {
    ht <- list(Herb=NULL, Tree=NULL)
    
    #plot x native & naturalized species matrix
    spe.nat <- natz.plots %>% filter(Habitat == q & INVASION.STATUS %in% c("native", "naturalized")) %>%
      mutate(presence = 1) %>%
      reshape2::dcast(., ID ~ Species, value.var="presence", fill = 0) %>% 
      column_to_rownames(var="ID")
    
    #beals co-occurrence scores
    b <- vegan::beals(spe.nat, include=F)*100
    bl <- SpecData.missFor[colnames(b), "Life.form"]
    bs <- SpecData.missFor[colnames(b), "INVASION.STATUS"]
    
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
      clusterExport(cl, list("null3ali.simulator", "apply", "as.randtest", "stats.ali", "stats.ali.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames", "Edistance_weighted"))
      clusterSetRNGStream(cl,1234)
      
      d <- pblapply(l, FUN = null3ali.simulator, cl = cl, 
                    beals.scores = b, bl = bl, bs = bs, tr = SpecData.missFor, tr.vars = tr.vars, nperm=999,
                    save.simulations = NULL)
      
      stopCluster(cl)
      
      ht[[m]] <- do.call(rbind.data.frame, d)
      ht[[m]] <- ht[[m]][ ht[[m]]$INVASION.STATUS != "native", ]
      
      rm(l)
    }
    
    NULL3ali.missFor[[q]] <- ht
  }
  
}

saveRDS(NULL3ali.missFor, file="NULL3ali.missFor.rds")