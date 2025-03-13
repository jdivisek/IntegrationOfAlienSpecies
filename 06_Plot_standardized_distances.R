#############################################
###    PLOT STANDARDIZED DISTANCES        ###
#############################################

#load packages
library(tidyverse)
library(introdataviz)
library(ggpubr)
library(rstatix)


###Standardized distances--------------------------------------------------

plot.names <- setNames(c("Grassland vegetation", "Ruderal and weed vegetation", 
                         "Rock and scree vegetation", "Wetland vegetation", 
                         "Scrub vegetation | Herbs & dwarf shrubs", "Scrub vegetation | Trees & tall shrubs",
                         "Forest vegetation | Herbs & dwarf shrubs", "Forest vegetation | Trees & tall shrubs"),
                       nm = c("T", "X", "S", "M", "KHerb", "KTree", "LHerb", "LTree"))

gg.plots <- list()

for(q in c("T", "X", "S", "M", "K", "L"))
{
  if(q %in% c("T", "X", "S", "M"))
  {
    dat1 <- STAT.NULLw[[q]] %>% filter(ID %in% unique(.$ID[.$INVASION.STATUS == "native"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "naturalized"])) %>% mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>% 
      filter(INVASION.STATUS %in% c("native", "naturalized")) %>% mutate(comparison = "natz.plots")

    dat2 <- STAT.NULLw[[q]] %>% filter(ID %in% unique(.$ID[.$INVASION.STATUS == "native"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "invasive"])) %>% mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>%
      filter(INVASION.STATUS %in% c("native", "invasive")) %>% mutate(comparison = "inv.plots")

    dat3 <- STATali.NULLw[[q]] %>% filter(ID %in% unique(.$ID[.$INVASION.STATUS == "naturalized"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "invasive"])) %>% mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>%
      filter(INVASION.STATUS %in% c("naturalized", "invasive")) %>% mutate(comparison = "ali.plots")

    dat_long <- rbind(dat1, dat2, dat3) %>% arrange(comparison, INVASION.STATUS, ID)
    dmax <- max(dat_long$`dm_O-E`)
    
    dat_long$INVASION.STATUS <- factor(dat_long$INVASION.STATUS, levels = c("native", "naturalized", "invasive"))
    dat_long$comparison <- factor(dat_long$comparison, levels = c("natz.plots", "inv.plots", "ali.plots"))
    
    ##Wilcoxon effect size and test
    eff <- dat_long %>% group_by(comparison) %>%  wilcox_effsize(data=., `dm_O-E` ~ INVASION.STATUS, paired = TRUE) %>% mutate(effsize = round(effsize, 2))
    test <- dat_long %>% group_by(comparison) %>% 
      rstatix::wilcox_test(data=., `dm_O-E` ~ INVASION.STATUS, p.adjust.method = "none", paired = T, detailed = F, alternative = "two.sided") %>% 
      adjust_pvalue(method = "fdr") %>% 
      add_significance("p.adj", cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")) %>% 
      mutate(p.adj = round(p.adj,3))

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
      scale_y_continuous(name = expression("Standardized distance (Δ"*italic(d)*")"),
                         expand = expansion(mult = c(0,0.18))) +
      scale_fill_manual(values = my.cols2) +
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
      annotate(geom="text", x = 1:3, y = dmax, label = paste("N =", test$n1), size=4.5, vjust = -1.7) +
      annotate(geom="text", x = 1:3, y = dmax, label = paste("r =", eff$effsize, test$p.adj.signif), size=4.5, vjust = 0.1) #+

  }
  
  if(q %in% c("K", "L") )
  {
    for(m in c("Herb", "Tree"))
    {
      dat1 <- STAT.NULLw[[q]][[m]] %>% filter(ID %in% unique(.$ID[.$INVASION.STATUS == "native"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "naturalized"])) %>% mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>% 
        filter(INVASION.STATUS %in% c("native", "naturalized")) %>% mutate(comparison = "natz.plots")

      dat2 <- STAT.NULLw[[q]][[m]] %>% filter(ID %in% unique(.$ID[.$INVASION.STATUS == "native"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "invasive"])) %>% mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>%
        filter(INVASION.STATUS %in% c("native", "invasive")) %>% mutate(comparison = "inv.plots")

      dat3 <- STATali.NULLw[[q]][[m]] %>% filter(ID %in% unique(.$ID[.$INVASION.STATUS == "naturalized"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "invasive"])) %>% mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>%
        filter(INVASION.STATUS %in% c("naturalized", "invasive")) %>% mutate(comparison = "ali.plots")

      dat_long <- rbind(dat1, dat2, dat3) %>% arrange(comparison, INVASION.STATUS, ID)
      dmax <- max(dat_long$`dm_O-E`)
      
      dat_long$INVASION.STATUS <- factor(dat_long$INVASION.STATUS, levels = c("native", "naturalized", "invasive"))
      dat_long$comparison <- factor(dat_long$comparison, levels = c("natz.plots", "inv.plots", "ali.plots"))
      
      ##Wilcoxon effect size and test
      eff <- dat_long %>% group_by(comparison) %>%  wilcox_effsize(data=., `dm_O-E` ~ INVASION.STATUS, paired = TRUE) %>% mutate(effsize = round(effsize, 2))
      test <- dat_long %>% group_by(comparison) %>% 
        rstatix::wilcox_test(data=., `dm_O-E` ~ INVASION.STATUS, p.adjust.method = "none", paired = T, detailed = F, alternative = "two.sided") %>% 
        adjust_pvalue(method = "fdr") %>% 
        add_significance("p.adj", cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")) %>% 
        mutate(p.adj = round(p.adj,3))

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
        scale_y_continuous(name = expression("Standardized distance (Δ"*italic(d)*")"),
                           expand = expansion(mult = c(0,0.18))) +
        scale_fill_manual(values = my.cols2) +
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
        annotate(geom="text", x = 1:3, y = dmax, label = paste("N =", test$n1), size=4.5, vjust = -1.7) +
        annotate(geom="text", x = 1:3, y = dmax, label = paste("r =", eff$effsize, test$p.adj.signif), size=4.5, vjust = 0.1) #+

    }
  }
}

tiff("Fig_2.tif", width = 11.5, height = 16, units = "in", compression = "lzw", res = 500)
ggarrange(plotlist = gg.plots, ncol=2, nrow=4)
dev.off()
