################################################
###       PRINCIPAL COMPONENT ANALYSIS       ###
################################################

library(factoextra)

##trait names
tr.vars <- c("HEIGHT", "SLA", "Germinule", "Leaf.area", "LDMC.mean", "FLOWERING.MEAN", "FLOWERING.LENGTH", "Gene.size.2C")

##habitat names
hab.names <- setNames(c("Grassland vegetation", "Ruderal and weed vegetation", "Rock and scree vegetation", "Wetland vegetation", "Scrub vegetation", "Forest vegetation"),
                      c("T", "X", "S", "M", "K", "L"))

##colors
my.cols2 <- c(native="#5EABB6",
              naturalized="goldenrod1",
              invasive="#A54888")

##PCA
for(q in names(hab.names))
{
  print(hab.names[q])
  
  sel <- spel.imp.ls %>% filter(Habitat == q) %>% select(Species) %>% .[,1] %>% unique()
  sel <- trait.imp.ls[rownames(trait.imp.ls) %in% sel, ]
  pca <- prcomp(sel[, tr.vars], center = F, scale. = F)
  rownames(pca$rotation) <- c("Plant height", "SLA", "Seed mass", "Leaf area",
                              "LDMC", "Middle of the\nflowering period", 
                              "Length of the\nflowering period", "Genome size")
  which(rownames(pca$x) == "Impatiens parviflora")
  
  if(q == "X"){
    pca$rotation[,1:2] <- pca$rotation[,1:2] *-1
    pca$x[,1:2] <- pca$x[,1:2] *-1}
  if(q == "S"){
    pca$rotation[,2] <- pca$rotation[,2] *-1
    pca$x[,2] <- pca$x[,2] *-1}
  if(q == "T"){
    pca$rotation[,1] <- pca$rotation[,1] *-1
    pca$x[,1] <- pca$x[,1] *-1}
  
  p1 <- fviz_pca_biplot(pca,
                        axes = c(1, 2),
                        geom.ind = "point",
                        geom.var = c("arrow", "text"),
                        repel = T,
                        palette = my.cols2,
                        habillage = factor(sel$INVASION.STATUS, levels=c("native", "naturalized", "invasive"), labels = c("native", "naturalized", "invasive")),
                        col.var = "black",
                        alpha.ind = 1,
                        addEllipses = T,
                        ellipse.level=0.95,
                        title = hab.names[q],
                        ggtheme = theme_bw(),
                        label.rectangle = T,
                        invisible = "quali",
                        labelsize = 3) + 
    theme(text = element_text(size = 16), legend.position = c(0.85,0.10), 
          legend.title = element_blank(), legend.text = element_text(size=14),
          plot.title = element_text(hjust = 0.5)) 
  
  pdf(paste0("PCA_", q, ".pdf"), width = 6.03, height = 5.92)
  plot(p1)
  dev.off()
}