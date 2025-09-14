###FIGURE 3---------------------------------------------------------------------
##Plot probability of overlap between native, naturalized and invasive species

##Load packages-----------------------------------------------------------------
library(ggpattern)
library(tidyverse)
library(introdataviz)
library(ggpubr)

##Define colors and names-------------------------------------------------------
my.cols <- setNames(c("#5EABB6", "goldenrod1", "#A54888"), nm = c("native", "naturalized", "invasive"))

plot.names <- setNames(c("Grassland vegetation", "Ruderal and weed vegetation", 
                         "Rock and scree vegetation", "Wetland vegetation", 
                         "Scrub vegetation | Herbs & dwarf shrubs", "Scrub vegetation | Trees & tall shrubs",
                         "Forest vegetation | Herbs & dwarf shrubs", "Forest vegetation | Trees & tall shrubs"),
                       nm = c("T", "X", "S", "M", "KHerb", "KTree", "LHerb", "LTree"))

##Plot Figure 3-----------------------------------------------------------------

gg.plots <- list()

for(q in c("T", "X", "S", "M", "K", "L"))
{
  if(q %in% c("T", "X", "S", "M"))
  {
    dat1 <- NULL3.missFor[[q]] %>% 
      filter(ID %in% unique(.$ID[.$INVASION.STATUS == "native"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "naturalized"])) %>% 
      mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>% 
      filter(INVASION.STATUS %in% c("naturalized")) %>% 
      mutate(comparison = "natz.plots")
    dat2 <- NULL3.missFor[[q]] %>% 
      filter(ID %in% unique(.$ID[.$INVASION.STATUS == "native"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "invasive"])) %>% 
      mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>%
      filter(INVASION.STATUS %in% c("invasive")) %>% 
      mutate(comparison = "inv.plots")
    dat3 <- NULL3ali.missFor[[q]] %>% 
      filter(ID %in% unique(.$ID[.$INVASION.STATUS == "naturalized"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "invasive"])) %>% 
      mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>%
      filter(INVASION.STATUS %in% c("invasive")) %>% 
      mutate(comparison = "ali.plots")
    
    dat_long <- rbind(dat1, dat2, dat3)
    head(dat_long)
    
    dat_long$INVASION.STATUS <- factor(dat_long$INVASION.STATUS, levels = c("native", "naturalized", "invasive"))
    dat_long$comparison <- factor(dat_long$comparison, levels = c("natz.plots", "inv.plots", "ali.plots"))

    ag <- aggregate(dat_long$`Edist_P-val`, by = list(dat_long$comparison), FUN = median)
    colnames(ag) <- c("comparison", "median")

    N <- table(dat_long$comparison)
    
    gg.plots[[q]] <- ggplot(dat_long, aes(x = comparison, y = `Edist_P-val`, fill = INVASION.STATUS)) +
      geom_flat_violin(position = position_nudge(x = 0), fill="gray", colour= NA, scale = "area",
                       trim=T, alpha = 0.6, show.legend = FALSE) +
      geom_boxplot_pattern(pattern = 'stripe', 
                           pattern_fill = my.cols[c(1,1,2)], 
                           pattern_colour = my.cols[c(2,2,3)], 
                           colour = "black",
                           pattern_spacing = 0.15,
                           pattern_density = 0.5,
                           pattern_linetype = 0, 
                           pattern_angle = 45,
                           pattern_alpha = 1,
                           position = position_nudge(x = -0.05), width = .3, alpha = 1, show.legend = FALSE, outliers = F) +
      stat_summary(fun="mean", geom = "point", show.legend = F, size = 4, shape=16, colour = "black",
                   position = position_nudge(x = -0.05)) +
      scale_x_discrete(name = "", labels = c("native vs. naturalized", "native vs. invasive", "naturalized vs. invasive")) +
      scale_y_continuous(name = "Probability",
                         breaks = seq(0, 1, 0.2), 
                         limits = c(0, 1.05)) +
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
      annotate(geom="text", x = 1, y = 1.05, label = paste("N =", N[[1]]), size=4.5) +
      annotate(geom="text", x = 2, y = 1.05, label = paste("N =", N[[2]]), size=4.5) +
      annotate(geom="text", x = 3, y = 1.05, label = paste("N =", N[[3]]), size=4.5)
    
  }
  if(q %in% c("K", "L") )
  {
    for(m in c("Herb", "Tree"))
    {
      dat1 <- NULL3.missFor[[q]][[m]] %>% 
        filter(ID %in% unique(.$ID[.$INVASION.STATUS == "native"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "naturalized"])) %>% 
        mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>% 
        filter(INVASION.STATUS %in% c("naturalized")) %>% 
        mutate(comparison = "natz.plots")
      dat2 <- NULL3.missFor[[q]][[m]] %>% 
        filter(ID %in% unique(.$ID[.$INVASION.STATUS == "native"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "invasive"])) %>% 
        mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>%
        filter(INVASION.STATUS %in% c("invasive")) %>% 
        mutate(comparison = "inv.plots")
      dat3 <- NULL3ali.missFor[[q]][[m]] %>% 
        filter(ID %in% unique(.$ID[.$INVASION.STATUS == "naturalized"]) & ID %in% unique(.$ID[.$INVASION.STATUS == "invasive"])) %>% 
        mutate(INVASION.STATUS = as.character(INVASION.STATUS)) %>%
        filter(INVASION.STATUS %in% c("invasive")) %>% 
        mutate(comparison = "ali.plots")

      dat_long <- rbind(dat1, dat2, dat3)
      head(dat_long)
      
      dat_long$INVASION.STATUS <- factor(dat_long$INVASION.STATUS, levels = c("native", "naturalized", "invasive"))
      dat_long$comparison <- factor(dat_long$comparison, levels = c("natz.plots", "inv.plots", "ali.plots"))
      
      ag <- aggregate(dat_long$`Edist_P-val`, by = list(dat_long$comparison), FUN = median)
      colnames(ag) <- c("comparison", "median")

      N <- table(dat_long$comparison)
      
      gg.plots[[paste0(q,m)]] <- ggplot(dat_long, aes(x = comparison, y = `Edist_P-val`, fill = INVASION.STATUS)) +
        geom_flat_violin(position = position_nudge(x = 0), fill="gray", colour= NA, scale = "area",
                         trim=T, alpha = 0.6, show.legend = FALSE) +
        geom_boxplot_pattern(pattern = 'stripe', 
                             pattern_fill = my.cols[c(1,1,2)], 
                             pattern_colour = my.cols[c(2,2,3)], 
                             colour = "black",
                             pattern_spacing = 0.15,
                             pattern_density = 0.5,
                             pattern_linetype = 0, 
                             pattern_angle = 45,
                             pattern_alpha = 1,
                             position = position_nudge(x = -0.05), width = .3, alpha = 1, show.legend = FALSE, outliers = F) +
        stat_summary(fun="mean", geom = "point", show.legend = F, size = 4, shape=16, colour="black",
                     position = position_nudge(x = -0.05)) +
        scale_x_discrete(name = "", labels = c("native vs. naturalized", "native vs. invasive", "naturalized vs. invasive")) +
        scale_y_continuous(name = "Probability",
                           breaks = seq(0, 1, 0.2), 
                           limits = c(0, 1.05)) +
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
        annotate(geom="text", x = 1, y = 1.05, label = paste("N =", N[[1]]), size=4.5) +
        annotate(geom="text", x = 2, y = 1.05, label = paste("N =", N[[2]]), size=4.5) +
        annotate(geom="text", x = 3, y = 1.05, label = paste("N =", N[[3]]), size=4.5)
      
    }
  }
}

tiff("Fig_3.tif", width = 11.5, height = 16, units = "in", compression = "lzw", res = 500)
# cairo_pdf("Fig_3.pdf", width = 11.5, height = 16)

ggarrange(plotlist = gg.plots, ncol=2, nrow=4)
dev.off()