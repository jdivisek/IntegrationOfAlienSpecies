##READ DATA---------------------------------------------------------------------

##Data for analyses available in the Zenodo repository:
#Divíšek, J., Pyšek, P., Richardson, D. M., Gotelli, N. J., Beckage, B., Molofsky, J., 
#Lososová, Z., & Chytrý, M. (2025). Data and R codes from: Naturalized and invasive 
#species integrate differently in the trait space of local plant communities (1.0) 
#[Data set]. Zenodo. https://doi.org/10.5281/zenodo.17116216

##Species data (invasion status, life form, frequency in habitat and traits)
#The file contains the following columns:
#
#NATIVE.ALIEN = Native/alien status of species (values: native, alien).
#INVASION.STATUS = Invasion status for alien species (values: native, naturalized, invasive).
#Life.form = Aggregated life forms.
#Herb.Tree = Herbs & dwarf shrubs vs. trees & tall shrubs.
#T = Species occurrence frequency (no. of vegetation plots) in grassland vegetation.
#X = Species occurrence frequency (no. of vegetation plots) in ruderal and weed vegetation.
#S = Species occurrence frequency (no. of vegetation plots) in rock and scree vegetation.
#M = Species occurrence frequency (no. of vegetation plots) in wetland vegetation.
#K = Species occurrence frequency (no. of vegetation plots) in scrub vegetation.
#L = Species occurrence frequency (no. of vegetation plots) in forest vegetation.
#HEIGHT = Maximum plant height (m).
#SLA = Specific leaf area (mm2 mg−1).
#Seed.mass = Seed mass (mg).
#Leaf.area = Leaf area (mm2).
#LDMC.mean = Leaf dry matter content (mg g-1).
#FLOWERING.MEAN = Middle of the flowering period (month).
#FLOWERING.LENGTH = Length of the flowering period (number of months). 
#Gen.size.2C = 2C genome size (Mbp).

#Dataset with missing trait values for some species
SpecData <- read.delim("SpecData.txt", header = T, row.names = 1)
head(SpecData)

#Missing trait values imputed by missForest function ("default" dataset)
SpecData.missFor <- read.delim("SpecData.missFor.txt", header = T, row.names = 1)

#Missing trait values imputed by MICE-PMM function
SpecData.MICE <- read.delim("SpecData.MICE.txt", header = T, row.names = 1)

#Missing trait values imputed by Phylopars function
SpecData.Rphylo <- read.delim("SpecData.Rphylo.txt", header = T, row.names = 1)

##Community data in long format
#The file contains the following columns:
#
#ID = ID of the vegetation plot.
#Species = Taxon name.
#cover = Taxon cover in the vegetation plot (%).
#Habitat = Habitat type (values: T, X, S, M, K, L)
#NATIVE.ALIEN = Native/alien status (values: native, alien).
#INVASION.STATUS = Invasion status for alien species (values: native, naturalized, invasive).
#Life.form = Aggregated life forms.
#Herb.Tree = Herbs & dwarf shrubs vs. trees & tall shrubs.

CommData <- read.delim("CommData.txt", header = T)
head(CommData)

##Prepare data for analyses----------------------------------------------------
library(tidyverse)
library(corrplot)

tr.vars <- c("HEIGHT", "SLA", "Seed.mass", "Leaf.area", "LDMC.mean", "FLOWERING.MEAN", "FLOWERING.LENGTH", "Gen.size.2C")

summary(SpecData.missFor[, tr.vars])
sapply(tr.vars, FUN = function(x) {hist(SpecData.missFor[,x], main =x)})

tr.cor <- round(cor(SpecData.missFor[, tr.vars], method = "spearman"), 3)

corrplot(tr.cor, method = 'color', diag = FALSE, type = 'lower', 
         tl.srt = 45, addCoef.col = 'black',
         pch.cex = 0.9, pch.col = 'grey20', tl.col = 'black')

CommData.missFor.raw <- cbind(CommData, SpecData.missFor[CommData$Species, tr.vars])

#transform and scale traits
SpecData.missFor[, tr.vars[-6]] <- log10(SpecData.missFor[, tr.vars[-6]])
SpecData.missFor[, tr.vars] <- scale(SpecData.missFor[, tr.vars])

CommData.missFor <- cbind(CommData, SpecData.missFor[CommData$Species, tr.vars])

##Plot trait spectrum for each habitat------------------------------------------

library(factoextra)

my.cols <- setNames(c("#5EABB6", "goldenrod1", "#A54888"), nm = c("native", "naturalized", "invasive"))
hab.names <- setNames(c("Grassland vegetation", "Ruderal and weed vegetation", 
                        "Rock and scree vegetation", "Wetland vegetation",
                        "Scrub vegetation", "Forest vegetation"),
                      nm = c("T", "X", "S", "M", "K", "L"))

for(q in names(hab.names))
{
  print(hab.names[q])
  
  sel <- SpecData.missFor[SpecData.missFor[, q] > 0, ]
  pca <- prcomp(sel[, tr.vars], center = F, scale. = F)
  rownames(pca$rotation) <- c("Plant height", "SLA", "Seed mass", "Leaf area",
                              "LDMC", "Middle of the\nflowering period", 
                              "Length of the\nflowering period", "Genome size")
  which(rownames(pca$x) == "Impatiens parviflora")
  
  if(q == "T"){
    pca$rotation[,1] <- pca$rotation[,1] *-1
    pca$x[,1] <- pca$x[,1] *-1}
  if(q == "S"){
    pca$rotation[,2] <- pca$rotation[,2] *-1
    pca$x[,2] <- pca$x[,2] *-1}

  p1 <- fviz_pca_biplot(pca,
                        axes = c(1, 2),
                        geom.ind = "point",
                        geom.var = c("arrow", "text"),
                        repel = T,
                        palette = my.cols,
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

##Middle of the flowering period as ordinal variable----------------------------

SpecData.missFor.ord <- read.delim("SpecData.missFor.txt", header = T, row.names = 1)

SpecData.missFor.ord[, tr.vars[-6]] <- log10(SpecData.missFor.ord[, tr.vars[-6]])
SpecData.missFor.ord[, tr.vars[-6]] <- scale(SpecData.missFor.ord[, tr.vars[-6]])
SpecData.missFor.ord$FLOWERING.MEAN <- ordered(SpecData.missFor.ord$FLOWERING.MEAN, levels = seq(3, 10, 0.5))


CommData.missFor.ord <- cbind(CommData, SpecData.missFor.ord[CommData$Species, tr.vars])
