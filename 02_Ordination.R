################################################
###       PRINCIPAL COMPONENT ANALYSIS       ###
################################################

library(factoextra)

##trait names
tr.vars <- c("HEIGHT", "SLA", "Germinule", "Leaf.area", "LDMC.mean", "FLOWERING.MEAN", "FLOWERING.LENGTH", "Gene.size.2C")

##colors
my.cols2 <- c(native="#5EABB6",
              naturalized="goldenrod1",
              invasive="#A54888")

##PCA
pca <- prcomp(trait.imp.ls[, tr.vars], center = T, scale. = T)

##PCA plot
p1 <- fviz_pca_biplot(pca,
                      axes = c(1, 2),
                      geom.ind = "point",
                      geom.var = c("arrow", "text"),
                      repel = T,
                      palette = my.cols2,
                      habillage = factor(trait.imp.ls$INVASION.STATUS, levels=c("native", "naturalized", "invasive"), labels = c("native", "naturalized", "invasive")),
                      col.var = "black",
                      alpha.ind = 0,
                      addEllipses = T,
                      title = "",
                      ggtheme = theme_bw(),
                      ellipse.type = "convex",
                      label.rectangle = T,
                      invisible = "quali",
                      labelsize = 3) + theme(text = element_text(size = 16), legend.position = c(0.85,0.15), legend.title = element_blank(), legend.text = element_text(size=16)) 
p1

##Find centroids
centroids <- aggregate(pca$x[,1:2], by=list(trait.imp.ls$INVASION.STATUS, trait.imp.ls$Life.form2), FUN=mean)
centroids$Group.1 <- factor(centroids$Group.1, levels=c("native", "naturalized", "invasive"), labels = c("native", "naturalized", "invasive"))
