##Simulates virtual communities by co-occurrence weighted drawing of species from the naturalized species pool.
#Beal's probability of species occurrence in the plot is used as a weight. The higher the probability, the higher the chance of being selected.
#This function uses Gower's (1971) distance with Podani's (1999) extension for
#ordinal variables. Gower distance matrix is then ordinated using principal
#coordinate analysis and species scores on ordination axes are used instead of 
#original trait values.
#Returns distance statistics, their deviations from null expectation and empirical p-values.

#@one.plot = Data frame with the species (rows) present in the community 
#and the following columns: ID, INVASION.STATUS, cover, Life.form, Herb.tree and traits.
#@beals.scores = Matrix of occurrence probabilities for the native and naturalized species in the vegetation plots returned by the vegan::beals function.
#@bl = Life forms of the native and naturalized species in the beals.scores matrix.
#@bs = Character vector of status labels (native and naturalized) for species in the beals.scores matrix.
#@tr = Data frame with species (rows) and their traits (columns).
#@tr.vars = Character vector of the trait names. Must match the column names in one.plot and tr.
#@nperm = Number of communities to be simulated.
#@save.simulations = Path to the folder in which the simulated communities are to be saved. The default value is NULL (simulated communities are not saved).

null3ali.gpcoa.simulator <- function(one.plot, beals.scores, bl, bs, tr, tr.vars, nperm=999, save.simulations=NULL)
{
  if(length(unique(one.plot$INVASION.STATUS)) > 2 & any(one.plot$INVASION.STATUS %in% c("native", "naturalized")))
  {
    op <- split(one.plot, one.plot$Life.form)##split based on life forms
    
    op.sim <- list()
    step <- 1
    repeat
    {
      op.sim[[step]] <- lapply(op, FUN = function(x)
      {
        ###sampling for alien species
        if(any(x$NATIVE.ALIEN == "alien"))
        {
          w <- setNames(beals.scores[as.character(x$ID[1]), bl==x$Life.form[1] & bs == "naturalized"],
                        colnames(beals.scores)[bl==x$Life.form[1] & bs == "naturalized"])
          if(sum(w) > 0){
            w <- w[w > 0] #remove species with Beals score < 0
          }
          else{ w[names(w)] <- 1 }
          
          size <- sum(x$INVASION.STATUS != "native")
          if(size <= length(w)) {sel <- sample(names(w), size = size, prob=w, replace=F)}
          else {sel <- sample(names(w), size = size, prob=w, replace=T)}
          
          x[x$INVASION.STATUS != "native", tr.vars] <- tr[sel, tr.vars]
          
          sel <- NULL
        }
        
        ###sampling for native species
        w <- setNames(beals.scores[as.character(x$ID[1]), bl==x$Life.form[1] & bs == "native"],
                      colnames(beals.scores)[bl==x$Life.form[1] & bs == "native"])
        if(sum(w) > 0){
          w <- w[w > 0] #remove species with Beals score < 0
        }
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
    
    obs <- stats.gpcoa.ali.calc(one.plot = one.plot, tr.vars = tr.vars)
    obs.t <- apply(X = one.plot[, tr.vars], MARGIN = 2, FUN = stats.gpcoa.ali, st = one.plot$INVASION.STATUS, 
                   cover = one.plot$cover, simplify = F)
    
    obs <- cbind(obs, do.call(cbind.data.frame, obs.t))
    
    ###list of simulations
    op.sim <- split(op.sim[,c("INVASION.STATUS", "cover", tr.vars)], op.sim$simulation)
    
    sim <- lapply(X = op.sim, FUN = function(x){
      si <- stats.gpcoa.ali(dat=x[, tr.vars], st=x$INVASION.STATUS, cover = x$cover)
      
      si.t <- apply(X = x[, tr.vars], MARGIN = 2, FUN = stats.gpcoa.ali, st = x$INVASION.STATUS, 
                    cover = x$cover, simplify = F)
      
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
