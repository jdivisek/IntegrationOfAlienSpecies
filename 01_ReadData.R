#############################
###       READ DATA       ###
#############################

##trait data with imputed values
trait.imp <- read.delim("trait_imp.txt", header=T, row.names=1, as.is = T)

#NATIVE.ALIEN = native & alien status
#INVASION.STATUS = invasion status
#HEIGHT = maximum height (m)
#SLA = specific lead area (mm2 mg−1)
#Germinule = seed weight (mg)
#Leaf.area = leaf area (mm2)
#LDMC.mean = leaf dry matter content (mg g-1)
#FLOWERING.MEAN = month in the middle of the flowering period (month)
#FLOWERING.LENGTH = length of the flowering period (months)
#Gene.size.2C = 2C genome size (Mbp)
#Life.form2 = broadly defined life forms
#Herb.Tree = Herbs & dwarf shrubs vs trees & tall shrubs
#T = species occurrence frequency in grassland vegetation
#X = species occurrence frequency in ruderal and weed vegetation
#S = species occurrence frequency in rock and scree vegetation
#M = species occurrence frequency in wetland vegetation
#K = species occurrence frequency in scrub vegetation
#L = species occurrence frequency in forest vegetation

##log-10 transformed and scaled trait values
trait.imp.ls <- read.delim("trait_imp_logScaled.txt", header=T, row.names=1, as.is = T)

##species composition data
spel.imp <- read.delim("species_imp.txt", header=T, row.names=1, as.is = T)

#ID = vegetation-plot ID
#Species = taxon name
#cover = taxon cover in %
#Habitat = plot classification to habitat types

##species composition data with log10-transformed and scaled traits
spel.imp.ls <- read.delim("species_imp_logScaled.txt", header=T, row.names=1, as.is = T)

###IDENTIFY PROMINENT INVADERS---------------------------------------------------

library(tidyverse)

out.spp <- list()

for(q in c("T", "X", "S", "M", "K", "L"))
{
  print(q)
  
  if(q %in% c("T", "X", "S", "M"))
  {
    ##naturalized
    sel <- spel.imp.ls.out %>% filter(Habitat == q & INVASION.STATUS == "naturalized") %>% 
      select(Species) %>% table() %>% sort(decreasing = T)
    sel <- sel[1:round(length(sel)*0.1,0)]
    out.spp[[q]][["naturalized"]] <- sel
    
    ##invasive
    sel <- spel.imp.ls.out %>% filter(Habitat == q & INVASION.STATUS == "invasive") %>% 
      select(Species) %>% table() %>% sort(decreasing = T)
    sel <- sel[1:round(length(sel)*0.1,0)]
    out.spp[[q]][["invasive"]] <- sel

  }
  if(q %in% c("K", "L"))
  {
    for(m in c("Herb", "Tree"))
    {
      ##naturalized
      sel <- spel.imp.ls.out %>% filter(Habitat == q & INVASION.STATUS == "naturalized" & Herb.Tree == m) %>%
        select(Species) %>% table() %>% sort(decreasing = T)
      sel <- sel[1:round(length(sel)*0.1,0)]
      out.spp[[q]][[m]][["naturalized"]] <- sel
      
      ##invasive
      sel <- spel.imp.ls.out %>% filter(Habitat == q & INVASION.STATUS == "invasive" & Herb.Tree == m) %>%
        select(Species) %>% table() %>% sort(decreasing = T)
      sel <- sel[1:round(length(sel)*0.1,0)]
      out.spp[[q]][[m]][["invasive"]] <- sel
      
    }
  }
  
}

out.spp

##PLOT PROMINENT INVADERS-------------------------------------------------------

tiff("Prominent_invaders.tif", 10, 28, units="in", res=500, compression = "lzw")

layout(matrix(1:16, ncol=2, byrow = T))
par(mar=c(9.5, 4, 2, 0.5), oma=rep(0.5, 4))

for(q in c("T", "X", "S", "M", "K", "L"))
{
  if(q %in% c("T", "X", "S", "M"))
  {
    ##naturalized
    sel <- spel.imp.ls %>% filter(Habitat == q & INVASION.STATUS == "naturalized") %>% 
      select(Species) %>% table() %>% sort(decreasing = T)
    sel <- sel[1:30]
    no <- length(out.spp[[q]][["naturalized"]])
    barplot(sel, col=c(rep(my.cols2[2],no), rep("gray", 30-no)), ylab = "No. of plots", las=2, cex.names=0.7, border = NA,
            main=paste(hab.names[q], "|", "naturalized spp."))
    
    ##invasive
    sel <- spel.imp.ls %>% filter(Habitat == q & INVASION.STATUS == "invasive") %>% 
      select(Species) %>% table() %>% sort(decreasing = T)
    sel <- sel[1:30]
    no <- length(out.spp[[q]][["invasive"]])
    barplot(sel, col=c(rep(my.cols2[3],no), rep("gray", 30-no)), ylab = "No. of plots", las=2, cex.names=0.7, border = NA,
            main=paste(hab.names[q], "|", "invasive spp."))
    
  }
  if(q %in% c("K", "L"))
  {
    for(m in c("Herb", "Tree"))
    {
      ##naturalized
      sel <- spel.imp.ls %>% filter(Habitat == q & INVASION.STATUS == "naturalized" & Herb.Tree == m) %>% 
        select(Species) %>% table() %>% sort(decreasing = T)
      sel <- sel[1:30]
      no <- length(out.spp[[q]][[m]][["naturalized"]])
      if(m == "Herb") {barplot(sel, col=c(rep(my.cols2[2],no), rep("gray", 30-no)), ylab = "No. of plots", las=2, cex.names=0.7, border = NA,
                               main=paste(hab.names[q], "| ", paste0(m, "s & dwarf shrubs"),  "|", "naturalized spp."))}
      if(m == "Tree") {barplot(sel, col=c(rep(my.cols2[2],no), rep("gray", 30-no)), ylab = "No. of plots", las=2, cex.names=0.7, border = NA,
                               main=paste(hab.names[q], "| ", paste0(m, "s & tall shrubs"),  "|", "naturalized spp."))}
      
      ##invasive
      sel <- spel.imp.ls %>% filter(Habitat == q & INVASION.STATUS == "invasive" & Herb.Tree == m) %>% 
        select(Species) %>% table() %>% sort(decreasing = T)
      sel <- sel[1:30]
      no <- length(out.spp[[q]][[m]][["invasive"]])
      if(m == "Herb"){barplot(sel, col=c(rep(my.cols2[3],no), rep("gray", 30-no)), ylab = "No. of plots", las=2, cex.names=0.7, border = NA,
                              main=paste(hab.names[q], "| ", paste0(m, "s & dwarf shrubs"),  "|", "invasive spp."))}
      if(m == "Tree"){barplot(sel, col=c(rep(my.cols2[3],no), rep("gray", 30-no)), ylab = "No. of plots", las=2, cex.names=0.7, border = NA,
                              main=paste(hab.names[q], "| ", paste0(m, "s & tall shrubs"),  "|", "invasive spp."))}
      
      
      
    }
  }
  
}

dev.off()



