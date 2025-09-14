##Wrapper for the stats.gpcoa.ali function. 
#Returns distance statistics for naturalized and invasive species in 
#the selected community (vegetation plot).
#This function uses Gower's (1971) distance with Podani's (1999) extension for
#ordinal variables. Gower distance matrix is then ordinated using principal
#coordinate analysis and species scores on ordination axes are used instead of 
#original trait values.

#@one.plot = Data frame with the species (rows) present in the community (vegetation plot) 
#and the following columns: ID, INVASION.STATUS, cover and traits.
#@tr.vars = Character vector of the trait names. Must match the column names in one.plot.

stats.gpcoa.ali.calc <- function(one.plot, tr.vars)
{
  if(length(unique(one.plot$INVASION.STATUS)) > 1 & any(one.plot$INVASION.STATUS %in% c("native", "naturalized")))
  {
    n <- as.data.frame(table(one.plot$INVASION.STATUS))
    colnames(n) <- c("INVASION.STATUS", "N")
    n <- cbind(n, stats.gpcoa.ali(dat = one.plot[,tr.vars], 
                                  st = one.plot$INVASION.STATUS,
                                  cover = one.plot$cover))
    n <- n %>% mutate(ID = one.plot$ID[1], .before=1)
    
    return(n)
  }
}
