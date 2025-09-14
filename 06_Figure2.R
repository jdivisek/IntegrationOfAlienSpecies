###FIGURE 2---------------------------------------------------------------------
##Plot mean distances (ΔD) of native, naturalized and invasive species from 
#the native center of the trait space in each plot

##Load packages-----------------------------------------------------------------
library(tidyverse)
library(introdataviz)
library(ggpubr)
library(rstatix)
library(coin)

##Define colors and names-------------------------------------------------------
my.cols <- setNames(c("#5EABB6", "goldenrod1", "#A54888"), nm = c("native", "naturalized", "invasive"))

plot.names <- setNames(c("Grassland vegetation", "Ruderal and weed vegetation", 
                         "Rock and scree vegetation", "Wetland vegetation", 
                         "Scrub vegetation | Herbs & dwarf shrubs", "Scrub vegetation | Trees & tall shrubs",
                         "Forest vegetation | Herbs & dwarf shrubs", "Forest vegetation | Trees & tall shrubs"),
                       nm = c("T", "X", "S", "M", "KHerb", "KTree", "LHerb", "LTree"))

##Plot Figure 2-----------------------------------------------------------------
gg.plots <- list()

for(q in c("T", "X", "S", "M", "K", "L"))
{
  if(q %in% c("T", "X", "S", "M"))
  {
    dat1 <- NULL3.missFor[[q]] %>% filter(ID %in% unique(.$ID[.$INVASION.STATUS == "native"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "naturalized"])) %>% 
      mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>% 
      filter(INVASION.STATUS %in% c("native", "naturalized")) %>% 
      mutate(comparison = "natz.plots")

    dat2 <- NULL3.missFor[[q]] %>% filter(ID %in% unique(.$ID[.$INVASION.STATUS == "native"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "invasive"])) %>% 
      mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>%
      filter(INVASION.STATUS %in% c("native", "invasive")) %>% 
      mutate(comparison = "inv.plots")

    dat3 <- NULL3ali.missFor[[q]] %>% filter(ID %in% unique(.$ID[.$INVASION.STATUS == "naturalized"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "invasive"])) %>% 
      mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>%
      filter(INVASION.STATUS %in% c("naturalized", "invasive")) %>% 
      mutate(comparison = "ali.plots")

    dat_long <- rbind(dat1, dat2, dat3) %>% arrange(comparison, INVASION.STATUS, ID)

    dmax <- max(dat_long$`dm_O-E`)
    
    dat_long$INVASION.STATUS <- factor(dat_long$INVASION.STATUS, levels = c("native", "naturalized", "invasive"))
    dat_long$comparison <- factor(dat_long$comparison, levels = c("natz.plots", "inv.plots", "ali.plots"))
    
    eff <- dat_long %>% group_by(comparison) %>%  
      wilcox_effsize(data=., `dm_O-E` ~ INVASION.STATUS, paired = TRUE) %>% 
      mutate(effsize = round(effsize, 2))
    test <- dat_long %>% group_by(comparison) %>% 
      rstatix::wilcox_test(data=., `dm_O-E` ~ INVASION.STATUS, p.adjust.method = "none", paired = T, detailed = F, alternative = "two.sided") %>% 
      adjust_pvalue(method = "fdr") %>% 
      add_significance("p.adj", cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")) %>% 
      mutate(p.adj = round(p.adj,3))

    ##Add z statistics
    eff[1, "z"] <- dat1 %>% 
      pivot_wider(id_cols = "ID",names_from = "INVASION.STATUS", values_from = "dm_O-E") %>% 
      coin::wilcoxsign_test(naturalized ~ native, data=., alternative = "two.sided", paired = T, zero.method="Pratt", distribution = "asymptotic") %>% 
      coin::statistic()
    eff[2, "z"] <- dat2 %>% 
      pivot_wider(id_cols = "ID",names_from = "INVASION.STATUS", values_from = "dm_O-E") %>% 
      coin::wilcoxsign_test(invasive ~ native, data=., alternative = "two.sided", paired = T, zero.method="Pratt", distribution = "asymptotic") %>% 
      coin::statistic()
    eff[3, "z"] <- dat3 %>% 
      pivot_wider(id_cols = "ID",names_from = "INVASION.STATUS", values_from = "dm_O-E") %>% 
      coin::wilcoxsign_test(invasive ~ naturalized, data=., alternative = "two.sided", paired = T, zero.method="Pratt", distribution = "asymptotic") %>% 
      coin::statistic()
    
    ##Add color for the direction of the effect
    eff$col <- sapply(eff$z, FUN = function(x){ifelse(x < 0, "red3", "dodgerblue3")})
    
    ##Significant results in bold
    test$fontface <- sapply(test$p.adj, FUN = function(x){ifelse(x < 0.05, 2, 1)})
    
    gg.plots[[q]] <- ggplot(dat_long, aes(x = comparison, y = `dm_O-E`, fill = INVASION.STATUS)) +
      geom_hline(yintercept=0, linetype="dashed", color = "red3") +
      introdataviz::geom_split_violin(alpha = .5, trim = T, colour= NA, scale="area", show.legend = F) +
      geom_boxplot(width = .25, alpha = 1, fatten = NULL, show.legend = FALSE,
                   position = position_dodge(.45), outliers = F) +
      stat_summary(fun="median", geom = "crossbar", show.legend = F, 
                   position = position_dodge(.45), width = 0.25) +
      stat_summary(fun="mean", geom = "point", show.legend = F, size = 3, shape=21,
                   position = position_dodge(.45)) +
      scale_x_discrete(name = "", labels = c("native vs. naturalized", "native vs. invasive", "naturalized vs. invasive")) +
      scale_y_continuous(name = "Deviation from null expectation (ΔD)",
                         expand = expansion(mult = c(0,0.27))) +
      scale_fill_manual(values = my.cols) +
      theme_bw() +
      theme(axis.text.y = element_text(size=12, angle = 90, hjust = 0.5, colour = "black"),
            axis.text.x = element_text(size=12, colour = "black"),
            axis.title.y = element_text(size=12),
            panel.grid.major.y = element_line(linetype = 3, size = 0.6, colour = "gray75"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
            plot.title = element_text(face = "bold", size = 15, hjust = 0.5)) +
      ggtitle(plot.names[q]) +
      annotate(geom="text", x = 1:3, y = dmax, label = paste("N =", test$n1), size=4.5, vjust = -3.3) +
      annotate(geom="text", x = 1:3, y = dmax, label = paste("z =", round(eff$z, 2), test$p.adj.signif), size=4.5, vjust = -1.6, colour = eff$col, fontface = test$fontface) +
      annotate(geom="text", x = 1:3, y = dmax, label = paste("r =", eff$effsize), size=4.5, vjust = 0.1)

  }
  if(q %in% c("K", "L"))
  {
    for(m in c("Herb", "Tree"))
    {
      dat1 <- NULL3.missFor[[q]][[m]] %>% 
        filter(ID %in% unique(.$ID[.$INVASION.STATUS == "native"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "naturalized"])) %>% 
        mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>% 
        filter(INVASION.STATUS %in% c("native", "naturalized")) %>% 
        mutate(comparison = "natz.plots")

      dat2 <- NULL3.missFor[[q]][[m]] %>% 
        filter(ID %in% unique(.$ID[.$INVASION.STATUS == "native"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "invasive"])) %>% 
        mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>%
        filter(INVASION.STATUS %in% c("native", "invasive")) %>% 
        mutate(comparison = "inv.plots")

      dat3 <- NULL3ali.missFor[[q]][[m]] %>% 
        filter(ID %in% unique(.$ID[.$INVASION.STATUS == "naturalized"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "invasive"])) %>% 
        mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>%
        filter(INVASION.STATUS %in% c("naturalized", "invasive")) %>% 
        mutate(comparison = "ali.plots")

      dat_long <- rbind(dat1, dat2, dat3) %>% arrange(comparison, INVASION.STATUS, ID)

      dmax <- max(dat_long$`dm_O-E`)
      
      dat_long$INVASION.STATUS <- factor(dat_long$INVASION.STATUS, levels = c("native", "naturalized", "invasive"))
      dat_long$comparison <- factor(dat_long$comparison, levels = c("natz.plots", "inv.plots", "ali.plots"))
      
      eff <- dat_long %>% 
        group_by(comparison) %>%  
        wilcox_effsize(data=., `dm_O-E` ~ INVASION.STATUS, paired = TRUE) %>% 
        mutate(effsize = round(effsize, 2))
      test <- dat_long %>% 
        group_by(comparison) %>% 
        rstatix::wilcox_test(data=., `dm_O-E` ~ INVASION.STATUS, p.adjust.method = "none", paired = T, detailed = F, alternative = "two.sided") %>% 
        adjust_pvalue(method = "fdr") %>% 
        add_significance("p.adj", cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")) %>% 
        mutate(p.adj = round(p.adj,3))

      ##Add z statistics
      eff[1, "z"] <- dat1 %>% 
        pivot_wider(id_cols = "ID",names_from = "INVASION.STATUS", values_from = "dm_O-E") %>% 
        coin::wilcoxsign_test(naturalized ~ native, data=., alternative = "two.sided", paired = T, zero.method="Pratt", distribution = "asymptotic") %>% 
        coin::statistic()
      eff[2, "z"] <- dat2 %>% 
        pivot_wider(id_cols = "ID",names_from = "INVASION.STATUS", values_from = "dm_O-E") %>% 
        coin::wilcoxsign_test(invasive ~ native, data=., alternative = "two.sided", paired = T, zero.method="Pratt", distribution = "asymptotic") %>% 
        coin::statistic()
      eff[3, "z"] <- dat3 %>% 
        pivot_wider(id_cols = "ID",names_from = "INVASION.STATUS", values_from = "dm_O-E") %>% 
        coin::wilcoxsign_test(invasive ~ naturalized, data=., alternative = "two.sided", paired = T, zero.method="Pratt", distribution = "asymptotic") %>% 
        coin::statistic()
      
      ##Add color for the direction of the effect
      eff$col <- sapply(eff$z, FUN = function(x){ifelse(x < 0, "red3", "dodgerblue3")})
      
      ##Significant results in bold
      test$fontface <- sapply(test$p.adj, FUN = function(x){ifelse(x < 0.05, 2, 1)})
      
      gg.plots[[paste0(q,m)]] <- ggplot(dat_long, aes(x = comparison, y = `dm_O-E`, fill = INVASION.STATUS)) +
        geom_hline(yintercept=0, linetype="dashed", color = "red3") +
        introdataviz::geom_split_violin(alpha = .5, trim = T, colour= NA, scale="area", show.legend = F) +
        geom_boxplot(width = .25, alpha = 1, fatten = NULL, show.legend = FALSE,
                     position = position_dodge(.45), outliers = F) +
        stat_summary(fun="median", geom = "crossbar", show.legend = F, 
                     position = position_dodge(.45), width = 0.25) +
        stat_summary(fun="mean", geom = "point", show.legend = F, size = 3, shape=21,
                     position = position_dodge(.45)) +
        scale_x_discrete(name = "", labels = c("native vs. naturalized", "native vs. invasive", "naturalized vs. invasive")) +
        scale_y_continuous(name = "Deviation from null expectation (ΔD)",
                           expand = expansion(mult = c(0,0.27))) +
        scale_fill_manual(values = my.cols) +
        theme_bw() +
        theme(axis.text.y = element_text(size=12, angle = 90, hjust = 0.5, colour = "black"),
              axis.text.x = element_text(size=12, colour = "black"),
              axis.title.y = element_text(size=12),
              panel.grid.major.y = element_line(linetype = 3, size = 0.6, colour = "gray75"),
              panel.grid.major.x = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
              plot.title = element_text(face = "bold", size = 15, hjust = 0.5)) +
        ggtitle(plot.names[paste0(q,m)]) +
        annotate(geom="text", x = 1:3, y = dmax, label = paste("N =", test$n1), size=4.5, vjust = -3.3) +
        annotate(geom="text", x = 1:3, y = dmax, label = paste("z =", round(eff$z, 2), test$p.adj.signif), size=4.5, vjust = -1.6, colour = eff$col, fontface = test$fontface) +
        annotate(geom="text", x = 1:3, y = dmax, label = paste("r =", eff$effsize), size=4.5, vjust = 0.1)
      
    }
  }
}

tiff("Fig_2.tif", width = 11.5, height = 16.5, units = "in", compression = "lzw", res = 500)
# cairo_pdf("Fig_2.pdf", width = 11.5, height = 16.5)

ggarrange(plotlist = gg.plots, ncol=2, nrow=4)
dev.off()