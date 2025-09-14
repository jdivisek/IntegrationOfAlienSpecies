##Simulates virtual communities by the unweighted drawing of species from the naturalized species pool.
#Returns distance statistics, their deviations from null expectation and empirical p-values 
#for alien species.

#@one.plot = Data frame with the species (rows) present in the community 
#and the following columns: ID, INVASION.STATUS, cover, Life.form, Herb.tree and traits.
#@pool.natv = Character vector of the native species present in the habitat species pool.
#@lf.natv = Character vector of the life forms of the native species present in the habitat species pool (must match the species names in pool.natv).
#@pool.natz = Character vector of the naturalized species present in the habitat species pool.
#@lf.natz = Character vector of the life forms of the naturalized species present in the habitat species pool (must match the species names in pool.natz).
#@tr = Data frame with species (rows) and their traits (columns).
#@tr.vars = Character vector of trait names. They must match the column names in one.plot and tr.
#@nperm = Number of communities to be simulated.
#@save.simulations = Path to the folder in which simulated communities are to be saved. The default value is NULL (simulated communities are not saved).

null1ali.simulator <- function(one.plot, pool.natv, lf.natv, pool.natz, lf.natz, tr, tr.vars, nperm=999, save.simulations=NULL)
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
          w <- pool.natz[lf.natz == x$Life.form[1]]
          
          size <- sum(x$NATIVE.ALIEN == "alien")
          if(size <= length(w)) {sel <- sample(w, size = size, replace=F)}
          else {sel <- sample(w, size = size, replace=T)}
          
          x[x$NATIVE.ALIEN == "alien", tr.vars] <- tr[sel, tr.vars]
          sel <- NULL
        }
        
        ###sampling for native species
        w <- pool.natv[lf.natv == x$Life.form[1]]
        
        size <- sum(x$NATIVE.ALIEN == "native")
        if(size <= length(w)) {sel <- sample(w, size = size, replace=F)}
        else {sel <- sample(w, size = size, replace=T)}
        
        x[x$NATIVE.ALIEN == "native", tr.vars] <- tr[sel, tr.vars]
        
        return(x)})
      
      op.sim[[step]] <- do.call(rbind.data.frame, op.sim[[step]]) %>% mutate(simulation = step)
      
      step <- step+1
      if(step == nperm+1)
      {
        break
      }
      
    }
    
    op.sim <- do.call(rbind.data.frame, op.sim)
    
    obs <- stats.ali.calc(one.plot = one.plot, tr.vars = tr.vars)
    obs.t <- apply(one.plot[, tr.vars], 2, FUN=function(x){
      stats.ali(x, st = one.plot$INVASION.STATUS, cover = one.plot$cover)}, simplify = F )
    obs <- cbind(obs, do.call(cbind.data.frame, obs.t))
    
    ###list of simulations
    op.sim <- split(op.sim[,c("INVASION.STATUS", "cover", tr.vars)], op.sim$simulation)
    
    sim <- lapply(op.sim, FUN = function(x){
      si <- stats.ali(dat=x[, tr.vars], st=x$INVASION.STATUS, cover = x$cover)
      
      si.t <- apply(x[, tr.vars], 2, FUN=function(y){
        stats.ali(y, st =x$INVASION.STATUS, cover = x$cover)}, simplify = F)
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
