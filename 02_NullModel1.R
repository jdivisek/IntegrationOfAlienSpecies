###NULL MODEL 1-----------------------------------------------------------------
##Simulates virtual communities by the unweighted drawing of species from the 
#native/naturalizes species pool.
#Returns distance statistics, their deviations from null expectation and 
#empirical p-values.

##Load packages-----------------------------------------------------------------
library(tidyverse)
library(ade4)
library(reshape2)
library(parallel)
library(pbapply)

##Load functions----------------------------------------------------------------

source("./Functions/Edistance_weighted.R")
source("./Functions/stats.R")
source("./Functions/stats.calc.R")
source("./Functions/stats.ali.R")
source("./Functions/stats.ali.calc.R")
source("./Functions/null1.simulator.R")
source("./Functions/null1ali.simulator.R")

##Run simulation for native vs alien (naturalized and invasive) species--------- 

#specify trait variables
tr.vars <- c("HEIGHT", "SLA", "Seed.mass", "Leaf.area", "LDMC.mean", "FLOWERING.MEAN", "FLOWERING.LENGTH", "Gen.size.2C")

NULL1.missFor <- list(T=list(), X=list(), S=list(), M=list(), K=list(), L=list())

for(q in c("T", "X", "S", "M", "K", "L"))
{
  print(q)

  if(q %in% c("T", "X", "S", "M"))
  {
    print("Running simulations...")

    l <- split(CommData.missFor[CommData.missFor$Habitat == q,], CommData.missFor[CommData.missFor$Habitat == q, "ID"])
    invaded <- unlist(lapply(l, FUN = function(x){length(unique(x$NATIVE.ALIEN))}))
    l <- l[invaded == 2]

    #sqrt transform covers
    l <- lapply(l, FUN = function(x) {
      x$cover <- sqrt(x$cover)
      return(x)})

    cl <- makeCluster(10) #number of cores for parallel computation
    clusterExport(cl, list("null1.simulator", "as.randtest", "stats", "stats.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames", "Edistance_weighted"))
    clusterSetRNGStream(cl,1234)

    d <- pblapply(l, FUN = null1.simulator, cl = cl,
                  pool = rownames(SpecData.missFor)[SpecData.missFor$NATIVE.ALIEN == "native" & SpecData.missFor[, q] > 0],
                  lf = SpecData.missFor$Life.form[SpecData.missFor$NATIVE.ALIEN == "native" & SpecData.missFor[, q] > 0],
                  tr = SpecData.missFor, tr.vars = tr.vars, nperm = 999,
                  save.simulations = NULL)

    stopCluster(cl)

    NULL1.missFor[[q]] <- do.call(rbind.data.frame, d)

    rm(l)
  }

  if(q %in% c("K", "L"))
  {
    ht <- list(Herb=NULL, Tree=NULL)
    
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

      cl <- makeCluster(10) #number of cores for parallel computation
      clusterExport(cl, list("null1.simulator", "as.randtest", "stats", "stats.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames", "Edistance_weighted"))
      clusterSetRNGStream(cl,1234)

      d <- pblapply(l, FUN = null1.simulator, cl = cl,
                    pool = rownames(SpecData.missFor)[SpecData.missFor$NATIVE.ALIEN == "native" & SpecData.missFor$Herb.Tree == m & SpecData.missFor[, q] > 0],
                    lf = SpecData.missFor$Life.form[SpecData.missFor$NATIVE.ALIEN == "native" & SpecData.missFor$Herb.Tree == m & SpecData.missFor[, q] > 0],
                    tr = SpecData.missFor, tr.vars = tr.vars, nperm = 999,
                    save.simulations = NULL)

      stopCluster(cl)

      ht[[m]] <- do.call(rbind.data.frame, d)

      rm(l)
    }

    NULL1.missFor[[q]] <- ht
  }

}

saveRDS(NULL1.missFor, file="NULL1.missFor.rds")

###Run simulation for naturalized vs invasive species--------------------------- 

NULL1ali.missFor <- list(T=list(), X=list(), S=list(), M=list(), K=list(), L=list())

#select plots with naturalized species
natz.id <- unique(CommData.missFor$ID[CommData.missFor$INVASION.STATUS == "naturalized"])
natz.plots <- CommData.missFor[CommData.missFor$ID %in% natz.id, ]

for(q in c("T", "X", "S", "M","K", "L"))
{
  print(q)

  if(q %in% c("T", "X", "S", "M"))
  {
    print("Running simulations...")

    l <- split(natz.plots[natz.plots$Habitat == q,], natz.plots[natz.plots$Habitat == q, "ID"])

    #sqrt transform covers
    l <- lapply(l, FUN = function(x) {
      x$cover <- sqrt(x$cover)
      return(x)})

    cl <- makeCluster(10)
    clusterExport(cl, list("null1ali.simulator", "as.randtest", "stats.ali", "stats.ali.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames", "Edistance_weighted"))
    clusterSetRNGStream(cl,1234)

    d <- pblapply(l, FUN = null1ali.simulator, cl = cl,
                  pool.natv = rownames(SpecData.missFor)[SpecData.missFor$NATIVE.ALIEN == "native" & SpecData.missFor[, q] > 0],
                  lf.natv = SpecData.missFor$Life.form[SpecData.missFor$NATIVE.ALIEN == "native" & SpecData.missFor[, q] > 0],
                  pool.natz = rownames(SpecData.missFor)[SpecData.missFor$INVASION.STATUS == "naturalized" & SpecData.missFor[, q] > 0],
                  lf.natz = SpecData.missFor$Life.form[SpecData.missFor$INVASION.STATUS == "naturalized" & SpecData.missFor[, q] > 0],
                  tr = SpecData.missFor, tr.vars = tr.vars, nperm = 999,
                  save.simulations = NULL)

    stopCluster(cl)

    NULL1ali.missFor[[q]] <- do.call(rbind.data.frame, d)
    NULL1ali.missFor[[q]] <- NULL1ali.missFor[[q]][ NULL1ali.missFor[[q]]$INVASION.STATUS != "native", ]

    rm(l)
  }

  if(q %in% c("K", "L"))
  {
    ht <- list(Herb=NULL, Tree=NULL)

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
      clusterExport(cl, list("null1ali.simulator", "as.randtest", "stats.ali", "stats.ali.calc", "%>%", "column_to_rownames", "mutate", "rownames_to_column", "remove_rownames", "Edistance_weighted"))
      clusterSetRNGStream(cl,1234)

      d <- pblapply(l, FUN = null1ali.simulator, cl = cl,
                    pool.natv = rownames(SpecData.missFor)[SpecData.missFor$NATIVE.ALIEN == "native" & SpecData.missFor$Herb.Tree == m & SpecData.missFor[, q] > 0],
                    lf.natv = SpecData.missFor$Life.form[SpecData.missFor$NATIVE.ALIEN == "native" & SpecData.missFor$Herb.Tree == m & SpecData.missFor[, q] > 0],
                    pool.natz = rownames(SpecData.missFor)[SpecData.missFor$INVASION.STATUS == "naturalized" & SpecData.missFor$Herb.Tree == m & SpecData.missFor[, q] > 0],
                    lf.natz = SpecData.missFor$Life.form[SpecData.missFor$INVASION.STATUS == "naturalized" & SpecData.missFor$Herb.Tree == m & SpecData.missFor[, q] > 0],
                    tr = SpecData.missFor, tr.vars = tr.vars, nperm = 999,
                    save.simulations = NULL)

      stopCluster(cl)

      ht[[m]] <- do.call(rbind.data.frame, d)
      ht[[m]] <- ht[[m]][ ht[[m]]$INVASION.STATUS != "native", ]

      rm(l)
    }

    NULL1ali.missFor[[q]] <- ht
  }

}

saveRDS(NULL1ali.missFor, file="NULL1ali.missFor.rds")