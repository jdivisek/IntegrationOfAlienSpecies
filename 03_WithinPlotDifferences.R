############################################################
###       WITHIN-PLOT TRAIT DIFFERENCES - HEAT MAP       ###
############################################################

##Calculate the differences between mean trait values of native and alien
##(either naturalized or invasive) species in each vegetation plot of each habitat
##type.
##Summarize within-plot differences for each habitat type using Cohen's d statistic

library(raster)
library(berryFunctions)
library(tidyverse)
library(reshape2)
library(misty)

##ASSEMBLE DATA-----------------------------------------------------------------
d <- list()
d$naturalized <- as.data.frame(matrix(data=NA, nrow=8, ncol=8, byrow = T))
colnames(d$naturalized) <- tr.vars
rownames(d$naturalized) <- c("T", "X", "S", "M", "KHerb", "KTree", "LHerb", "LTree")

d$invasive <- as.data.frame(matrix(data=NA, nrow=8, ncol=8, byrow = T))
colnames(d$invasive) <- tr.vars
rownames(d$invasive) <- c("T", "X", "S", "M", "KHerb", "KTree", "LHerb", "LTree")

for(q in c("T", "X", "S", "M", "K", "L"))
{
  if(q %in% c("T", "X", "S", "M"))
  {
    sel <- spel.imp.ls %>% dplyr::filter(Habitat == q) %>% select(all_of(c("ID", "INVASION.STATUS", tr.vars)))
    ag <- aggregate(sel[, tr.vars], by=list(sel$ID, sel$INVASION.STATUS), FUN = mean)
    colnames(ag)[1:2] <- c("ID", "INVASION.STATUS")
    
    for(tr in tr.vars)
    {
      tab <- dcast(ag, ID ~ INVASION.STATUS, value.var = tr, fill = NA)
      
      d$naturalized[q, tr] <- cohens.d(tab$native, tab$naturalized, paired = T, weighted = F, cor=T, correct = T, na.omit=T)$result$d
      d$invasive[q, tr] <- cohens.d(tab$native, tab$invasive, paired = T, weighted = F, cor=T, correct = T, na.omit=T)$result$d
    }
    
  }
  if(q %in% c("K", "L"))
  {
    for(m in c("Herb", "Tree"))
    {
      sel <- spel.imp.ls %>% dplyr::filter(Habitat == q & Herb.Tree == m) %>% select(all_of(c("ID", "INVASION.STATUS", tr.vars)))
      ag <- aggregate(sel[, tr.vars], by=list(sel$ID, sel$INVASION.STATUS), FUN = mean)
      colnames(ag)[1:2] <- c("ID", "INVASION.STATUS")
      
      for(tr in tr.vars)
      {
        tab <- dcast(ag, ID ~ INVASION.STATUS, value.var = tr, fill = NA)
        
        d$naturalized[paste0(q, m), tr] <- cohens.d(tab$native, tab$naturalized, paired = T, weighted = F, cor=T, correct = T, na.omit=T)$result$d
        d$invasive[paste0(q, m), tr] <- cohens.d(tab$native, tab$invasive, paired = T, weighted = F, cor=T, correct = T, na.omit=T)$result$d
      }
      
    }
    
  }
}

rownames(d$naturalized) <- c(hab.names[1:4], paste(hab.names[c(5,5,6,6)], "|", c("Herbs & dwarf shrubs", "Trees & tall shrubs")))
rownames(d$invasive) <- c(hab.names[1:4], paste(hab.names[c(5,5,6,6)], "|", c("Herbs & dwarf shrubs", "Trees & tall shrubs")))
colnames(d$naturalized) <- c("Plant height", "SLA", "Seed weight", "Leaf area",  "LDMC", "Middle of the flowering period", "Length of the flowering period", "Genome size")
colnames(d$invasive) <- c("Plant height", "SLA", "Seed weight", "Leaf area",  "LDMC", "Middle of the flowering period", "Length of the flowering period", "Genome size")

##PLOT HEATMAP------------------------------------------------------------------

plot.data <- list()

r <- raster(as.matrix(d$naturalized), xmn=0, xmx=8, ymn=0, ymx=8)
plot.data$naturalized <- rasterToPoints(r)

r <- raster(as.matrix(d$invasive), xmn=0, xmx=8, ymn=0, ymx=8)
plot.data$invasive <- rasterToPoints(r)

tiff("Heatmap_Cohens.d.tif", width = 13.5, height = 7, units = "cm", res = 400, compression = "lzw")

m <- matrix(data=0, ncol=(8*2)+1, nrow=8+2)
m[1:8, 1:17] <- 1
m[9:10, 1:17] <- 2

par(mar=c(0,0,0,0), oma=c(0,8,7,1), xaxs="i", yaxs="i")
layout(m)
# layout.show()

range(c(plot.data$naturalized[,3], plot.data$invasive[,3]))

plot(plot.data$naturalized[,1:2], type='n', axes=F, xlab="", ylab="", asp=1, ylim=c(0,8), xlim=c(0,17))
br <- c(-3.5, -0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8, 3.5)

my.ramp.cold <- colorRampPalette(c("white", "#1A6179")); my.ramp.hot <- colorRampPalette(c("white", "#F88C24"))
cols <- c(rev(my.ramp.cold(5)[-c(1)]), my.ramp.hot(5)[-c(1)])

#naturalized species
clas <- cut(plot.data$naturalized[,3], br, include.lowest = TRUE, labels=1:(length(br)-1))
colPoints(x=plot.data$naturalized[,1], y=plot.data$naturalized[,2], z=plot.data$naturalized[,3], col=cols,
          method = "custom", breaks=br, pch=15, cex=3.4, axes=F, add=T, legend = F, asp=1)
text(x=rep(0, 8), y=seq(0.5,7.5), labels=rev(rownames(d$naturalized)), cex=0.7, xpd=NA, pos=2)
text(x=seq(0.5, 7.5), y=rep(8.2,8), labels=colnames(d$naturalized) , cex=0.7, xpd=NA, adj=c(0,0), srt=45)
text(4, 12.5,labels="Naturalized species",cex=1, xpd=NA, font=2)

#invasive species
clas <- cut(plot.data$invasive[,3], br, include.lowest = TRUE, labels=1:(length(br)-1))
colPoints(x=plot.data$invasive[,1]+9, y=plot.data$invasive[,2], z=plot.data$invasive[,3], col=cols,
          method = "custom", breaks=br ,pch=15, cex=3.4, add=T, legend = F, asp=1)
text(x=seq(0.5, 7.5)+9, y=rep(8.2,8), labels= colnames(d$invasive), cex=0.7, xpd=NA, adj=c(0,0), srt=45)
text(4+9, 12.5,labels="Invasive species",cex=1, xpd=NA, font=2)

#plot legend
plot(1, type="n", axes=F, xlab="", ylab="")

colPointsLegend(z=c(-4, -3, -2, -1, 0, 1, 2, 3, 4),
                col=cols, x1=0.2, x2=0.8, y1=0.4, y2=0.65, density = F, lines=F, 
                title=expression("Cohen's"~italic(d)~"for the within-plot differences between native and alien species"), cex=0.6, atminmax=F, resetfocus = F, mar=c(0,0,0,0),
                bb = c(-4, -3, -2, -1, 0, 1, 2, 3, 4), 
                at = c(-3, -2, -1, 0, 1, 2, 3), 
                labels = c(-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8))

dev.off()

###BOXPLOTS---------------------------------------------------------------------

##Boxplots of the within-plot trait differences

##replace d value by "*" or "+" signs
d.sign <- function(v, sign="*"){
  v <- abs(v)
  s <- rep("", times=length(v))
  s[v >= 0.2] <- sign
  s[v >= 0.5] <- paste0(rep(sign,2), collapse = "")
  s[v >= 0.8] <- paste0(rep(sign,3), collapse = "")
  return(s)
}

tr.names <- c("Plant height", "SLA", "Seed weight", "Leaf area",
              "LDMC", "Middle of the flowering period", "Length of the flowering period",
              "Genome size")
hab.names <- setNames(c("Grassland vegetation", "Ruderal and weed vegetation",
                        "Rock and scree vegetation", "Wetland vegetation",
                        "Scrub vegetation", "Forest vegetation"), c("T", "X", "S", "M", "K", "L"))

tiff("Within-plot_differences.tif", 6, 7, units="in", res=500, compression = "lzw")
layout(matrix(1:8, ncol=2, byrow = T))
par(mar=c(0.5,0.5,1.5,0.5), oma=c(5,2.5,0,0), mgp=c(2, 0.6, 0))

for(q in c("T", "X", "S", "M", "K", "L"))
{
  print(q)
  if(q %in% c("T", "X", "S", "M"))
  {
    plot(NA, xlim= c(0.5,8.5), ylim=c(-3,4), xaxt = "n", yaxt = "n", xlab="", ylab="",
         main=hab.names[q], cex.main=0.8, cex.axis=0.6, cex.lab=0.8)
    
    if(q %in% c("T", "S")){
      axis(2, cex.axis=0.8)
      text(-0.7,0.5, labels=expression("Δ"*bar(italic("T"))), srt=90, xpd=NA, cex=0.8)}
    
    
    sel <- spel.imp.ls %>% filter(Habitat == q) %>% select(all_of(c("ID", "INVASION.STATUS", tr.vars)))
    sel <- aggregate(sel[,tr.vars], by=list(sel$ID, sel$INVASION.STATUS), mean)
    
    for(i in 1:length(tr.vars))
    {
      tab <- dcast(sel, Group.1 ~ Group.2, value.var = tr.vars[i], fill = NA)[,-1]
      
      boxplot(tab$naturalized - tab$native, border="black", col=adjustcolor(my.cols2[2], 1),
              at = i-0.2, boxwex=0.6, axes = F, add=T, outline=F, staplewex=0, whisklty="solid", xlab="", ylab="")
      points(i-0.2, mean(tab$naturalized - tab$native, na.rm=T), pch=21, bg=my.cols2[2], cex=1)
      
      boxplot(tab$invasive - tab$native, border="black", col=adjustcolor(my.cols2[3], 1),
              at = i+0.2, boxwex=0.6, axes = F, add=T, outline=F, staplewex=0, whisklty="solid", xlab="", ylab="")
      points(i+0.2, mean(tab$invasive - tab$native, na.rm=T), pch=21, bg=my.cols2[3], cex=1)
      
      
      d <- cohens.d(tab$native, tab$naturalized, paired = T, weighted = F, cor=T, correct = F, na.omit=T)$result$d
      text(x=i, y=3.2, labels = ifelse(d > 0, d.sign(d, sign="+"), d.sign(d, sign="-")), cex=1.1)
      d <- cohens.d(tab$native, tab$invasive, paired = T, weighted = F, cor=T, correct = T, na.omit=T)$result$d
      text(x=i, y=3.7, labels = ifelse(d > 0, d.sign(d, sign="+"), d.sign(d, sign="-")), cex=1.1)
      
    }
    
    abline(h=0, col=my.cols2[1], lwd=1, lty="dashed")
    
  }
  if(q %in% c("K", "L"))
  {
    for(m in c("Herb", "Tree"))
    {
      if(m == "Herb"){plot(NA, xlim= c(0.5,8.5), ylim=c(-3,4), xaxt = "n", yaxt = "n", xlab="", ylab="",
                           main=paste(hab.names[q], "|", paste0(m, "s & dwarf shrubs")), cex.main=0.8, cex.axis=0.6, cex.lab=0.8)}
      if(m == "Tree"){plot(NA, xlim= c(0.5,8.5), ylim=c(-3,4), xaxt = "n", yaxt = "n", xlab="", ylab="",
                           main=paste(hab.names[q], "|", paste0(m, "s & tall shrubs")), cex.main=0.8, cex.axis=0.6, cex.lab=0.8)}
      
      if(m == "Herb"){
        axis(2, cex.axis=0.8)
        text(-0.7,0.5, labels=expression("Δ"*bar(italic("T"))), srt=90, xpd=NA, cex=0.8)}
      
      sel <- spel.imp.ls %>% filter(Habitat == q & Herb.Tree == m) %>% select(all_of(c("ID", "INVASION.STATUS", tr.vars)))
      sel <- aggregate(sel[,tr.vars], by=list(sel$ID, sel$INVASION.STATUS), mean)
      
      for(i in 1:length(tr.vars))
      {
        tab <- dcast(sel, Group.1 ~ Group.2, value.var = tr.vars[i], fill = NA)[,-1]
        
        boxplot(tab$naturalized - tab$native, border="black", col=adjustcolor(my.cols2[2], 1),
                at = i-0.2, boxwex=0.6, axes = F, add=T, outline=F, staplewex=0, whisklty="solid", xlab="", ylab="")
        points(i-0.2, mean(tab$naturalized - tab$native, na.rm=T), pch=21, bg=my.cols2[2], cex=1)
        
        boxplot(tab$invasive - tab$native, border="black", col=adjustcolor(my.cols2[3], 1),
                at = i+0.2, boxwex=0.6, axes = F, add=T, outline=F, staplewex=0, whisklty="solid", xlab="", ylab="")
        points(i+0.2, mean(tab$invasive - tab$native, na.rm=T), pch=21, bg=my.cols2[3], cex=1)
        
        
        d <- cohens.d(tab$native, tab$naturalized, paired = T, weighted = F, cor=T, correct = T, na.omit=T)$result$d
        text(x=i, y=3.2, labels = ifelse(d > 0, d.sign(d, sign="+"), d.sign(d, sign="-")), cex=1.1)
        d <- cohens.d(tab$native, tab$invasive, paired = T, weighted = F, cor=T, correct = T, na.omit=T)$result$d
        text(x=i, y=3.7, labels = ifelse(d > 0, d.sign(d, sign="+"), d.sign(d, sign="-")), cex=1.1)
        
      }
      
      
      if(q == "L"){
        axis(side = 1, at=1:8, labels = FALSE)
        text(x = 1:8, y = par("usr")[3] - 0.90,
             labels = tr.names, xpd = NA,
             srt = 25, adj = 0.965, cex = 0.8)}
      abline(h=0, col=my.cols2[1], lwd=1, lty="dashed")
      
    }
  }
}

dev.off()

