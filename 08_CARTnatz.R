##REGRESSION TREE FOR NATURALIZED SPECIES---------------------------------------

##Load packages-----------------------------------------------------------------
library(rpart)
library(maptree)
library(partykit)
library(tidyverse)

##Assemble data-----------------------------------------------------------------
sel <- list(NULL3.missFor[["T"]], NULL3.missFor[["X"]], NULL3.missFor[["S"]], NULL3.missFor[["M"]], NULL3.missFor[["K"]][["Herb"]],  NULL3.missFor[["L"]][["Herb"]]) %>% 
  do.call(rbind.data.frame, .) %>% 
  dplyr::select(all_of(c("ID", "N", "INVASION.STATUS", "dm", "dm_O-E"))) %>% 
  mutate(ID2 = paste0(ID, INVASION.STATUS))

sel <- CommData.missFor.raw %>% 
  filter(Herb.Tree == "Herb") %>% 
  dplyr::select(all_of(c("ID", "Habitat"))) %>% 
  left_join(sel, ., by="ID", multiple = "first")

ag <- CommData.missFor.raw %>% 
  filter(Herb.Tree == "Herb") %>% 
  dplyr::select(all_of(c("ID", "INVASION.STATUS", tr.vars))) %>% 
  group_by(ID, INVASION.STATUS) %>% 
  summarise(HEIGHT= mean(HEIGHT), SLA=mean(SLA), Seed.mass=mean(Seed.mass), 
            Leaf.area= mean(Leaf.area), LDMC.mean=mean(LDMC.mean), 
            FLOWERING.MEAN = mean(FLOWERING.MEAN), FLOWERING.LENGTH = mean(FLOWERING.LENGTH), 
            Gen.size.2C= mean(Gen.size.2C)) %>% 
  ungroup() %>% 
  mutate(ID2 = paste0(ID, INVASION.STATUS)) %>% 
  dplyr::select(all_of(c("ID2", tr.vars))) %>% 
  as.data.frame()

cart.data <- left_join(sel, ag, by="ID2") %>% filter(INVASION.STATUS == "naturalized")

cart.data$`dm_O-E` <- round(cart.data$`dm_O-E`, 3)
cart.data$HEIGHT <- round(cart.data$HEIGHT, 1)
cart.data$SLA <- round(cart.data$SLA, 0)
cart.data$Seed.mass <- round(cart.data$Seed.mass, 1)
cart.data$Leaf.area <- round(cart.data$Leaf.area, 0)
cart.data$LDMC.mean <- round(cart.data$LDMC.mean, 0)
cart.data$FLOWERING.MEAN <- round(cart.data$FLOWERING.MEAN, 1)
cart.data$FLOWERING.LENGTH <- round(cart.data$FLOWERING.LENGTH, 1)
cart.data$Gen.size.2C <- round(cart.data$Gen.size.2C,0)

##CART--------------------------------------------------------------------------
set.seed(1234)
tree1.natz <- rpart(`dm_O-E` ~ HEIGHT + SLA + Seed.mass + Leaf.area + LDMC.mean + FLOWERING.MEAN + FLOWERING.LENGTH + Gen.size.2C, 
               data = cart.data, xval=10) 
summary(tree1.natz)
plotcp(tree1.natz)

windows()
draw.tree(tree1.natz, nodeinfo = TRUE, digits = 1, cex =0.8)

tree2.natz <- as.party(tree1.natz) 
plot(tree2.natz, cex=0.8)