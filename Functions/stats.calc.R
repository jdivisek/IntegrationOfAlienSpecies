##Wrapper for the stats function. 
#Returns distance statistics between native and alien (either naturalized or invasive) species in 
#the selected community (vegetation plot)

#@one.plot = Data frame with the species (rows) present in the community (vegetation plot) 
#and the following columns: ID, INVASION.STATUS, cover and traits.
#@tr.vars = Character vector of the trait names. Must match the column names in one.plot.

stats.calc <- function(one.plot, tr.vars)
{
  if(length(unique(one.plot$INVASION.STATUS)) > 1 & any(one.plot$INVASION.STATUS %in% c("native")))
  {
    n <- as.data.frame(table(one.plot$INVASION.STATUS))
    colnames(n) <- c("INVASION.STATUS", "N")
    n <- cbind(n, stats(dat = one.plot[,tr.vars], 
                        st = one.plot$INVASION.STATUS,
                        cover = one.plot$cover))
    n <- n %>% mutate(ID = one.plot$ID[1], .before=1)
    
    return(n)
  }
}