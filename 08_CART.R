##################################
###       DECISION TREES       ###
##################################

#load packages
library(rpart)
library(maptree)
library(partykit)
library(tidyverse)

##Select herbs and dwarf shrubs
sel <- list(ACC.NULL3.2[["T"]], ACC.NULL3.2[["X"]], ACC.NULL3.2[["S"]], ACC.NULL3.2[["M"]], ACC.NULL3.2[["K"]][["Herb"]],  ACC.NULL3.2[["L"]][["Herb"]]) %>% 
  do.call(rbind.data.frame, .) %>% select(all_of(c("ID", "N", "INVASION.STATUS", "acc", "acc.obs.minus.exp"))) %>% 
  mutate(ID2 = paste0(ID, INVASION.STATUS))
head(sel)

sel <- spel.imp.ip %>% filter(Herb.Tree == "Herb") %>% select(all_of(c("ID", "Habitat"))) %>% left_join(sel, ., by="ID", multiple = "first")

ag <- spel.imp.ip %>% filter(Herb.Tree == "Herb") %>% select(all_of(c("ID", "INVASION.STATUS", tr.vars))) %>% 
  group_by(ID, INVASION.STATUS) %>% summarise(HEIGHT= mean(HEIGHT), SLA=mean(SLA), Germinule=mean(Germinule),
                                              Leaf.area= mean(Leaf.area), LDMC.mean=mean(LDMC.mean), FLOWERING.MEAN = mean(FLOWERING.MEAN),
                                              FLOWERING.LENGTH = mean(FLOWERING.LENGTH), Gene.size.2C= mean(Gene.size.2C)) %>% ungroup() %>% 
  mutate(ID2 = paste0(ID, INVASION.STATUS)) %>% select(all_of(c("ID2", tr.vars))) %>% as.data.frame()

dat1 <- left_join(sel, ag, by="ID2")
head(dat1)

#Select trees and tall shrubs
sel <- list(ACC.NULL3.2[["K"]][["Tree"]],  ACC.NULL3.2[["L"]][["Tree"]]) %>% 
  do.call(rbind.data.frame, .) %>% select(all_of(c("ID", "N", "INVASION.STATUS", "acc", "acc.obs.minus.exp"))) %>% 
  mutate(ID2 = paste0(ID, INVASION.STATUS))
head(sel)

sel <- spel.imp %>% filter(Herb.Tree == "Tree") %>% select(all_of(c("ID", "Habitat"))) %>% left_join(sel, ., by="ID", multiple = "first")

ag <- spel.imp %>% filter(Herb.Tree == "Tree") %>% select(all_of(c("ID", "INVASION.STATUS", tr.vars))) %>% 
  group_by(ID, INVASION.STATUS) %>% summarise(HEIGHT= mean(HEIGHT), SLA=mean(SLA), Germinule=mean(Germinule),
                                              Leaf.area= mean(Leaf.area), LDMC.mean=mean(LDMC.mean), FLOWERING.MEAN = mean(FLOWERING.MEAN),
                                              FLOWERING.LENGTH = mean(FLOWERING.LENGTH), Gene.size.2C= mean(Gene.size.2C)) %>% ungroup() %>% 
  mutate(ID2 = paste0(ID, INVASION.STATUS)) %>% select(all_of(c("ID2", tr.vars))) %>% as.data.frame()

dat2 <- left_join(sel, ag, by="ID2")
head(dat2)

##with trees
plot.data <- rbind(dat1[dat1$INVASION.STATUS != "native", ],
                   dat2[dat2$INVASION.STATUS != "native", ]) 

##no trees
plot.data <- dat1[dat1$INVASION.STATUS != "native", ]

##trees only
plot.data <- dat2[dat2$INVASION.STATUS != "native", ]

###CART
set.seed(1234)
tree1 <- rpart(acc.obs.minus.exp ~ HEIGHT + SLA + Germinule + Leaf.area + LDMC.mean + FLOWERING.MEAN + FLOWERING.LENGTH + Gene.size.2C, 
               data = plot.data, xval=10) 
summary(tree1)
plotcp(tree1)

windows()
draw.tree(tree1, nodeinfo = TRUE, digits = 1, cex =0.8)

tree2 <- as.party(tree1) 
plot(tree2)

