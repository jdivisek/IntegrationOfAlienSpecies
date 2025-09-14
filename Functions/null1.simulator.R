##Simulates virtual communities by the unweighted drawing of species from the native species pool.
#Returns distance statistics, their deviations from null expectation and empirical p-values.

#@one.plot = Data frame with the species (rows) present in the community 
#and the following columns: ID, INVASION.STATUS, cover, Life.form, Herb.tree and traits.
#@pool = Character vector of the native species present in the habitat species pool.
#@lf = Character vector of the life forms of native species present in the habitat species pool (must match the species names in the pool).
#@tr = Data frame with species (rows) and their traits (columns).
#@tr.vars = Character vector of trait names. They must match the column names in one.plot and tr.
#@nperm = Number of communities to be simulated.
#@save.simulations = Path to the folder in which simulated communities are to be saved. The default value is NULL (simulated communities are not saved).

null1.simulator <- function(one.plot, pool, lf, tr, tr.vars, nperm=999, save.simulations=NULL)
{
  if(length(unique(one.plot$INVASION.STATUS)) > 1 & any(one.plot$INVASION.STATUS %in% c("native")))
  {
    op <- split(one.plot, one.plot$Life.form)##split based on life forms
    
    op.sim <- list()
    step <- 1
    repeat
    {
      op.sim[[step]] <- lapply(op, FUN = function(x)
      {
        w <- pool[lf == x$Life.form[1]]
        
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
    
    obs <- stats.calc(one.plot = one.plot, tr.vars = tr.vars)
    obs.t <- apply(X = one.plot[, tr.vars], MARGIN = 2, FUN = function(x){stats(x, st = one.plot$INVASION.STATUS, cover = one.plot$cover)}, simplify = FALSE)
    obs <- cbind(obs, do.call(cbind.data.frame, obs.t))
    
    ###list of simulations
    op.sim <- split(op.sim[,c("INVASION.STATUS", "cover", tr.vars)], op.sim$simulation)
    
    sim <- lapply(X = op.sim, FUN = function(x){
      si <- stats(dat=x[, tr.vars], st=x$INVASION.STATUS, cover = x$cover)
      
      si.t <- apply(X = x[, tr.vars], MARGIN = 2, FUN = function(y){stats(y, st =x$INVASION.STATUS, cover = x$cover)}, simplify = FALSE)
      si <- cbind(si, do.call(cbind.data.frame, si.t))
      
      si <- mutate(.data = si, INVASION.STATUS = rownames(si), .before = 1)
      
      return(si[order(rownames(si)), ])
    })
    sim <- do.call(rbind.data.frame, sim)
    
    merged <- rbind(obs[, !colnames(obs) %in% c("ID", "N")], sim)
    
    for(st in obs$INVASION.STATUS)
    {
      ## Observed minus expected values
      obs[obs$INVASION.STATUS == st, paste0(colnames(merged)[-1], "_O-E")] <- apply(X = merged[merged$INVASION.STATUS == st, -1], MARGIN = 2, FUN = function(x){ x[1] - mean(x[-1])}) 
      
      ##Empirical p-value for distances
      obs[obs$INVASION.STATUS == st, paste0(colnames(merged)[-1], "_P-val")] <- apply(X = merged[merged$INVASION.STATUS == st, -1], MARGIN = 2, FUN = function(x){ 
        as.randtest(sim = x[-1], obs = x[1], alter = "greater")$pvalue}) 
      
    }
    
    if(!is.null(save.simulations))
    {
      write.table(sim, file=paste0(save.simulations, "/", one.plot$ID[1], ".txt"), sep="\t", dec=".")
    }
    
    return(obs)
    
  }
  
}
