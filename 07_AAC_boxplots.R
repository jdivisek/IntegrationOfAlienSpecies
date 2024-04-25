#############################################
###        BOXPLOTS FOR ACC VALUES        ###
#############################################

#load packages
library(misty)
library(tidyverse)
library(reshape2)

###MULTIDIMENSIONAL AAC VALUES--------------------------------------------------

hab.names <- setNames(c("Grassland vegetation",
                        "Ruderal and weed vegetation",
                        "Rock and scree vegetation",
                        "Wetland vegetation",
                        "Scrub vegetation",
                        "Forest vegetation"), c("T", "X", "S", "M", "K", "L"))

tiff("BoxACC.NULLw.tif", width =  11.5, height =  6.5, units = "in", res = 500, compression = "lzw")

par(mar= c(2, 2, 2.5, 0.5), oma=c(0.5,2,1,0))
layout(matrix(1:8, nrow=2, ncol = 4, byrow = T))

for(q in c("T", "X", "S", "M", "K", "L"))#
{
  if(q %in% c("T", "X", "S", "M"))
  {
    ses <- dcast(ACC.NULLw[[q]], ID ~ INVASION.STATUS, value.var="acc.obs.minus.exp", fill = NA) %>% column_to_rownames(var="ID") %>% .[, c("native", "naturalized", "invasive")]
    
    ymax <- max(boxplot(ses, data=ses, plot=F)$stats[5,])
    ymin <- min(boxplot(ses, data=ses, plot=F)$stats[1,])
    if(q %in% c("T", "X")) {ymax <- ymin+((ymax-ymin)*1.3)}
    else {ymax <- ymin+((ymax-ymin)*1.2)}
    
    boxplot(ses, xaxt="n", yaxt="n", border=F, axes=F, outline=F, ylim=c(ymin, ymax),
            main=hab.names[q], col=NA)
    grid()
    box()
    
    # jittered points
    x.point <- seq(0.62, 1.38, length.out = length(ses$native))
    points(sample(x.point), ses$native, pch=16, cex=0.3, col=my.cols2[1])

    x.point <- seq(1.62, 2.38, length.out = length(ses$naturalized))
    points(sample(x.point), ses$naturalized, pch=16, cex=0.3, col=my.cols2[2])

    x.point <- seq(2.62, 3.38, length.out = length(ses$invasive))
    points(sample(x.point), ses$invasive, pch=16, cex=0.3, col=my.cols2[3])

    boxplot(ses, outline=FALSE, add=T, axes=F,
            col = c(adjustcolor(my.cols2[1], 0.5),
                    adjustcolor(my.cols2[2], 0.5),
                    adjustcolor(my.cols2[3], 0.5)), 
            lwd=1, notch=F, ylim=c(ymin, ymax), staplewex=0, whisklty="solid")
    
    ##plot means 
    points(1, mean(ses$native, na.rm=T), pch=21, cex=2, bg=my.cols2[1])
    points(2, mean(ses$naturalized, na.rm=T), pch=21, cex=2, bg=my.cols2[2])
    points(3, mean(ses$invasive, na.rm=T), pch=21, cex=2, bg=my.cols2[3])
    
    abline(h=0, col="red3", lty="solid", lwd=2)
    # abline(h=c(-1.96,1.96), col="black", lty="dashed")
    axis(1, at=1:3, labels = c("native", "naturalized", "invasive"))
    axis(2)
    if(q == "T"){
      text(-0.1, par("usr")[3]+(diff(par("usr")[3:4])/2), label="observed AAC - expected AAC", xpd=NA, srt=90)
    }
    
    text(1, ymax, paste0("N = ", sum(!is.na(ses$native))) ,pos=1, cex=1.1)
    text(2, ymax, paste0("N = ", sum(!is.na(ses$naturalized))) ,pos=1, cex=1.1)
    text(3, ymax, paste0("N = ", sum(!is.na(ses$invasive))) ,pos=1, cex=1.1)
    
    ###add Cohen's d
    d <- cohens.d(ses$native, ses$naturalized, paired = T, weighted = F, cor=T, correct = T, na.omit=T)$result$d
    d <- round(d, 2)
    text(2, ymin+((ymax-ymin)*0.94), bquote(italic("d")~"="~.(d)), pos=1, cex=1.1)
    d <- cohens.d(ses$native, ses$invasive, paired = T, weighted = F, cor=T, correct = T, na.omit=T)$result$d
    d <- round(d, 2)
    text(3, ymin+((ymax-ymin)*0.94), bquote(italic("d")~"="~.(d)), pos=1, cex=1.1)
    
  }
  
  if(q %in% c("K", "L"))
  {
    for(m in c("Herb", "Tree"))
    {
      ses <- dcast(ACC.NULLw[[q]][[m]], ID ~ INVASION.STATUS, value.var="acc.obs.minus.exp", fill = NA) %>% column_to_rownames(var="ID") %>% .[, c("native", "naturalized", "invasive")]
      
      ymax <- max(boxplot(ses, data=ses, plot=F)$stats[5,])
      ymin <- min(boxplot(ses, data=ses, plot=F)$stats[1,])
      if(q %in% c("T", "X")) {ymax <- ymin+((ymax-ymin)*1.3)}
      else {ymax <- ymin+((ymax-ymin)*1.2)}
      
      if(m == "Herb"){boxplot(ses, xaxt="n", yaxt="n", border=F, axes=F, outline=F, ylim=c(ymin, ymax),
                              main=paste(hab.names[q], "|", paste0(m, "s & dwarf shrubs")), col=NA)}
      if(m == "Tree"){boxplot(ses, xaxt="n", yaxt="n", border=F, axes=F, outline=F, ylim=c(ymin, ymax),
                              main=paste(hab.names[q], "|", paste0(m, "s & tall shrubs")), col=NA)}
      grid()
      box()
      
      # jittered points
      x.point <- seq(0.62, 1.38, length.out = length(ses$native))
      points(sample(x.point), ses$native, pch=16, cex=0.3, col=my.cols2[1])

      x.point <- seq(1.62, 2.38, length.out = length(ses$naturalized))
      points(sample(x.point), ses$naturalized, pch=16, cex=0.3, col=my.cols2[2])

      x.point <- seq(2.62, 3.38, length.out = length(ses$invasive))
      points(sample(x.point), ses$invasive, pch=16, cex=0.3, col=my.cols2[3])

      boxplot(ses, outline=FALSE, add=T, axes=F,
              col = c(adjustcolor(my.cols2[1], 0.5),
                      adjustcolor(my.cols2[2], 0.5),
                      adjustcolor(my.cols2[3], 0.5)), 
              lwd=1, notch=F, ylim=c(ymin, ymax), staplewex=0, whisklty="solid")
      
      ##plot means 
      points(1, mean(ses$native, na.rm=T), pch=21, cex=2, bg=my.cols2[1])
      points(2, mean(ses$naturalized, na.rm=T), pch=21, cex=2, bg=my.cols2[2])
      points(3, mean(ses$invasive, na.rm=T), pch=21, cex=2, bg=my.cols2[3])
      
      abline(h=0, col="red3", lty="solid", lwd=2)
      axis(1, at=1:3, labels = c("native", "naturalized", "invasive"))
      axis(2)
      if(q == "K" & m == "Herb"){
        text(-0.1, par("usr")[3]+(diff(par("usr")[3:4])/2), label="observed AAC - expected AAC", xpd=NA, srt=90)
      }
      
      text(1, ymax, paste0("N = ", sum(!is.na(ses$native))) ,pos=1, cex=1.1)
      text(2, ymax, paste0("N = ", sum(!is.na(ses$naturalized))) ,pos=1, cex=1.1)
      text(3, ymax, paste0("N = ", sum(!is.na(ses$invasive))) ,pos=1, cex=1.1)
      
      ###add Cohen's d
      d <- cohens.d(ses$native, ses$naturalized, paired = T, weighted = F, cor=T, correct = T, na.omit=T)$result$d
      d <- round(d, 2)
      text(2, ymin+((ymax-ymin)*0.94), bquote(italic("d")~"="~.(d)), pos=1, cex=1.1)
      d <- cohens.d(ses$native, ses$invasive, paired = T, weighted = F, cor=T, correct = T, na.omit=T)$result$d
      d <- round(d, 2)
      text(3, ymin+((ymax-ymin)*0.94), bquote(italic("d")~"="~.(d)), pos=1, cex=1.1)
      
    }
  }
}

dev.off()

###AAC VALUES FOR INDIVIDUAL TRAITS---------------------------------------------

tr.names <- c("Plant height", "SLA", "Seed weight", "Leaf area",
              "LDMC", "Middle of the flowering period", "Length of the flowering period",
              "Genome size")
hab.names <- setNames(c("Grassland vegetation", "Ruderal and weed vegetation",
                        "Rock and scree vegetation", "Wetland vegetation",
                        "Scrub vegetation", "Forest vegetation"), c("T", "X", "S", "M", "K", "L"))

tiff("BoxACC.trait.NULLw.tif", 6, 7, units="in", res=500, compression = "lzw")
layout(matrix(1:8, ncol=2, byrow = T))
par(mar=c(0.5,0.5,1.5,0.5), oma=c(5,2.5,0,0), mgp=c(2, 0.6, 0))

for(q in c("T", "X", "S", "M", "K", "L"))
{
  print(q)
  if(q %in% c("T", "X", "S", "M"))
  {
    plot(NA, xlim= c(0.5,8.5), ylim=c(-0.6,0.9), xaxt = "n", yaxt = "n", xlab="", ylab="",
         main=hab.names[q], cex.main=0.8, cex.axis=0.6, cex.lab=0.8)
    polygon(c(par("usr")[1], 1.5, 1.5, rev(par("usr")[1])),
            c(par("usr")[c(2,2)], par("usr")[c(3,3)]), col="gray95", border=NA)
    polygon(c(2.5, 3.5, 3.5, 2.5),
            c(par("usr")[c(2,2)], par("usr")[c(3,3)]), col="gray95", border=NA)
    polygon(c(4.5, 5.5, 5.5, 4.5),
            c(par("usr")[c(2,2)], par("usr")[c(3,3)]), col="gray95", border=NA)
    polygon(c(6.5, 7.5, 7.5, 6.5),
            c(par("usr")[c(2,2)], par("usr")[c(3,3)]), col="gray95", border=NA)
    box()
    
    if(q %in% c("T", "S")){
      axis(2, cex.axis=0.8)
      text(-0.7, mean(par("usr")[3:4]), labels="observed AAC - expected AAC", srt=90, xpd=NA, cex=0.8)}
    
    abline(h=0, col="red")
    
    sel <- ACC.NULLw[[q]] %>% filter(INVASION.STATUS == "native") %>% select(all_of(c("ID", "INVASION.STATUS", paste0(tr.vars, ".obs.minus.exp")))) %>% 
      gather(trait, value, HEIGHT.obs.minus.exp:Gene.size.2C.obs.minus.exp, factor_key=TRUE) 
    boxplot(value ~ trait, data=sel, border="black", col=adjustcolor(my.cols2[1], 1),
            at = c(1:8)-0.25, boxwex=0.15, axes = F, add=T, outline=F, staplewex=0, whisklty="solid", xlab="", ylab="", lwd=0.8)
    points(c(1:8)-0.25, aggregate(sel$value, by=list(sel$trait), mean)$x, pch=21, bg=my.cols2[1], cex=0.8, lwd=0.8)
    
    sel <- ACC.NULLw[[q]] %>% filter(INVASION.STATUS == "naturalized") %>% select(all_of(c("ID", "INVASION.STATUS", paste0(tr.vars, ".obs.minus.exp")))) %>% 
      gather(trait, value, HEIGHT.obs.minus.exp:Gene.size.2C.obs.minus.exp, factor_key=TRUE)
    boxplot(value ~ trait, data=sel, border="black", col=adjustcolor(my.cols2[2], 1),
            at = c(1:8), boxwex=0.15, axes = F, add=T, outline=F, staplewex=0, whisklty="solid", xlab="", ylab="", lwd=0.8)
    points(c(1:8), aggregate(sel$value, by=list(sel$trait), mean)$x, pch=21, bg=my.cols2[2], cex=0.8, lwd=0.8)
    
    sel <- ACC.NULLw[[q]] %>% filter(INVASION.STATUS == "invasive") %>% select(all_of(c("ID", "INVASION.STATUS", paste0(tr.vars, ".obs.minus.exp")))) %>% 
      gather(trait, value, HEIGHT.obs.minus.exp:Gene.size.2C.obs.minus.exp, factor_key=TRUE)
    boxplot(value ~ trait, data=sel, border="black", col=adjustcolor(my.cols2[3], 1),
            at = c(1:8)+0.25, boxwex=0.15, axes = F, add=T, outline=F, staplewex=0, whisklty="solid", xlab="", ylab="", lwd=0.8)
    points(c(1:8)+0.25, aggregate(sel$value, by=list(sel$trait), mean)$x, pch=21, bg=my.cols2[3], cex=0.8, lwd=0.8)
    
    ##get Cohen's d
    d.natz <- list(); d.inv <- list()
    for(tr in tr.vars)
    {
      tab <- dcast(ACC.NULLw[[q]], ID ~ INVASION.STATUS, value.var = paste0(tr, ".obs.minus.exp"), fill=NA)
      d.natz[[tr]] <- cohens.d(tab$native, tab$naturalized, paired = T, weighted = F, cor=T, correct = T, na.omit=T)$result$d
      d.inv[[tr]] <- cohens.d(tab$native, tab$invasive, paired = T, weighted = F, cor=T, correct = T, na.omit=T)$result$d
    }
    
    text(x=1:8, y=0.77, labels = unlist(lapply(d.natz, FUN = function(x){ifelse(x > 0, d.sign(x, sign="+"), d.sign(x, sign="-"))})), cex=1.1)
    text(x=1:8, y=0.87, labels = unlist(lapply(d.inv, FUN = function(x){ifelse(x > 0, d.sign(x, sign="+"), d.sign(x, sign="-"))})), cex=1.1)
    
  }
  if(q %in% c("K", "L"))
  {
    for(m in c("Herb", "Tree"))
    {
      if(m == "Herb"){plot(NA, xlim= c(0.5,8.5), ylim=c(-0.6,0.9), xaxt = "n", yaxt = "n", xlab="", ylab="",
                           main=paste(hab.names[q], "|", paste0(m, "s & dwarf shrubs")), cex.main=0.8, cex.axis=0.6, cex.lab=0.8)}
      if(m == "Tree"){plot(NA, xlim= c(0.5,8.5), ylim=c(-0.6,0.9), xaxt = "n", yaxt = "n", xlab="", ylab="",
                           main=paste(hab.names[q], "|", paste0(m, "s & tall shrubs")), cex.main=0.8, cex.axis=0.6, cex.lab=0.8)}
      
      polygon(c(par("usr")[1], 1.5, 1.5, rev(par("usr")[1])),
              c(par("usr")[c(2,2)], par("usr")[c(3,3)]), col="gray95", border=NA)
      polygon(c(2.5, 3.5, 3.5, 2.5),
              c(par("usr")[c(2,2)], par("usr")[c(3,3)]), col="gray95", border=NA)
      polygon(c(4.5, 5.5, 5.5, 4.5),
              c(par("usr")[c(2,2)], par("usr")[c(3,3)]), col="gray95", border=NA)
      polygon(c(6.5, 7.5, 7.5, 6.5),
              c(par("usr")[c(2,2)], par("usr")[c(3,3)]), col="gray95", border=NA)
      box()
      
      
      if(m == "Herb"){
        axis(2, cex.axis=0.8)
        text(-0.7,mean(par("usr")[3:4]), labels="observed AAC - expected AAC", srt=90, xpd=NA, cex=0.8)}
      abline(h=0, col="red")
      
      sel <- ACC.NULLw[[q]][[m]] %>% filter(INVASION.STATUS == "native") %>% select(all_of(c("ID", "INVASION.STATUS", paste0(tr.vars, ".obs.minus.exp")))) %>% 
        gather(trait, value, HEIGHT.obs.minus.exp:Gene.size.2C.obs.minus.exp, factor_key=TRUE) 
      boxplot(value ~ trait, data=sel, border="black", col=adjustcolor(my.cols2[1], 1),
              at = c(1:8)-0.25, boxwex=0.15, axes = F, add=T, outline=F, staplewex=0, whisklty="solid", xlab="", ylab="", lwd=0.8)
      points(c(1:8)-0.25, aggregate(sel$value, by=list(sel$trait), mean)$x, pch=21, bg=my.cols2[1], cex=0.8, lwd=0.8)
      
      sel <- ACC.NULLw[[q]][[m]] %>% filter(INVASION.STATUS == "naturalized") %>% select(all_of(c("ID", "INVASION.STATUS", paste0(tr.vars, ".obs.minus.exp")))) %>% 
        gather(trait, value, HEIGHT.obs.minus.exp:Gene.size.2C.obs.minus.exp, factor_key=TRUE)
      boxplot(value ~ trait, data=sel, border="black", col=adjustcolor(my.cols2[2], 1),
              at = c(1:8), boxwex=0.15, axes = F, add=T, outline=F, staplewex=0, whisklty="solid", xlab="", ylab="", lwd=0.8)
      points(c(1:8), aggregate(sel$value, by=list(sel$trait), mean)$x, pch=21, bg=my.cols2[2], cex=0.8, lwd=0.8)
      
      sel <- ACC.NULLw[[q]][[m]] %>% filter(INVASION.STATUS == "invasive") %>% select(all_of(c("ID", "INVASION.STATUS", paste0(tr.vars, ".obs.minus.exp")))) %>% 
        gather(trait, value, HEIGHT.obs.minus.exp:Gene.size.2C.obs.minus.exp, factor_key=TRUE)
      boxplot(value ~ trait, data=sel, border="black", col=adjustcolor(my.cols2[3], 1),
              at = c(1:8)+0.25, boxwex=0.15, axes = F, add=T, outline=F, staplewex=0, whisklty="solid", xlab="", ylab="", lwd=0.8)
      points(c(1:8)+0.25, aggregate(sel$value, by=list(sel$trait), mean)$x, pch=21, bg=my.cols2[3], cex=0.8, lwd=0.8)
      
      if(q == "L"){
        axis(side = 1, at=1:8, labels = FALSE)
        text(x = 1:8, y = par("usr")[3] - 0.18,
             labels = tr.names, xpd = NA,
             srt = 25, adj = 0.965, cex = 0.8)}
      
      ##get Cohen's d
      d.natz <- list(); d.inv <- list()
      for(tr in tr.vars)
      {
        tab <- dcast(ACC.NULLw[[q]][[m]], ID ~ INVASION.STATUS, value.var = paste0(tr, ".obs.minus.exp"), fill=NA)
        d.natz[[tr]] <- cohens.d(tab$native, tab$naturalized, paired = T, weighted = F, cor=T, correct = T, na.omit=T)$result$d
        d.inv[[tr]] <- cohens.d(tab$native, tab$invasive, paired = T, weighted = F, cor=T, correct = T, na.omit=T)$result$d
      }
      text(x=1:8, y=0.77, labels = unlist(lapply(d.natz, FUN = function(x){ifelse(x > 0, d.sign(x, sign="+"), d.sign(x, sign="-"))})), cex=1.1)
      text(x=1:8, y=0.87, labels = unlist(lapply(d.inv, FUN = function(x){ifelse(x > 0, d.sign(x, sign="+"), d.sign(x, sign="-"))})), cex=1.1)
      
    }
  }
}

dev.off()

