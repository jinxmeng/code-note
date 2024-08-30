#### Info ####
# Encoding: utf-8
# Author: Jinxin Meng
# Email: mengjx855@163.com
# Created Data：2023-9-15
# Modified Data: 2023-9-21
# Version: 1.0

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

#### calcu_PCA ####
# PCA analysis
# input: profile | distance
# dim: how many dimension is used to analyze
# cumulative_eig: 轴累积解释度达到多少阈值，输出多少个轴
# prefix: 修改列名的前缀
# add_eig: 列名中含有解释度
calcu_PCA <- function(profile, dim = 2, cumulative_eig = NULL, prefix = NULL, add_eig = F) {
  profile <- data.frame(profile, check.names = F)
  PCA <- prcomp(t(profile), scale. = T, center = T)
  PCA_sum <- summary(PCA)
  PCA_eig <- round(PCA_sum$importance[2,]*100, 2)
  if (is.null(cumulative_eig)) {
    PCA_points <- data.frame(PCA$x[,(1:dim)]) %>% 
      rownames_to_column(var = "sample")
  } else if (!is.null(cumulative_eig)) {
    tmp_vec <- cumsum(PCA_eig)
    for (i in 1:length(tmp_vec)) {
      if (tmp_vec[i] >= cumulative_eig) {
        dim <- i
        break
      }
    }
    PCA_points <- data.frame(PCA$x[,(1:dim)]) %>%
      rownames_to_column(var = "sample")
  }
  
  if (!is.null(prefix)) {
    prefix <- as.character(prefix)
    colnames(PCA_points)[2:(dim + 1)] <- paste0(prefix, "_", seq(1, dim))
  } else {
    colnames(PCA_points)[2:(dim + 1)] <- paste0("PCA", seq(1, dim))
  }
  
  if (isTRUE(add_eig)) {
    colnames(PCA_points)[2:(dim + 1)] <- paste0(colnames(PCA_points)[2:(dim + 1)], "(", PCA_eig[1:dim], "%)")
  }
  
  out = list(points = PCA_points, 
             dim = dim, 
             eig = PCA_eig, 
             eig_ = paste0(colnames(PCA_points)[2:(dim + 1)], "(", PCA_eig[1:dim], "%)"))
  
  return(out)
}

#### plot_PCA ####
# input --> exp table & sample_group
# group_order --> group level
# group_color --> group color
# file format
# exp|profile|otu
#         s1  s2 ...
# gene1   2   3  
# gene2   12  2  
# ...
# group
# sampple group        
# gene1   ctr
# gene2   case
# ...
plot_PCA <- function(exp, group, group_order = NULL, group_color = NULL, 
                     display_type = "line") {
  if(is.null(group_order)){ group_order <- unique(group$group) } else { group_order <- group_order }
  if(is.null(group_color)){
    color <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3",
               "#a6d854","#ffd92f","#e5c494","#b3b3b3")
    group_color <- rep(color, times = ceiling(length(group_order)/length(color)))[seq_len(length(group_order))]
  } else { group_color <- group_color }
  # PCA
  PCA <- prcomp(t(exp), scale. = T, center = T)
  PCA_points <- data.frame(PCA$x[,1:2]) %>% 
    rownames_to_column(var = "sample")
  PCA_sum <- summary(PCA)
  xlab <- paste0("PC1 (", round(PCA_sum$importance[2,1]*100, 2), "%)")
  ylab <- paste0("PC2 (", round(PCA_sum$importance[2,2]*100, 2), "%)")
  plotdat <- merge(PCA_points, group, by = "sample", all.x = T) %>% 
    mutate(group = factor(group, levels = group_order))
  if (display_type == "line") {
    plotdat_mean <- plotdat %>% 
      dplyr::select(-sample) %>% 
      group_by(group) %>% 
      summarise_all(mean) %>% 
      dplyr::rename(X1mean = PC1, X2mean = PC2) %>% 
      merge(x = plotdat %>% 
              dplyr::select(-sample) %>% 
              relocate(group), 
            y = ., by = "group") %>% 
      mutate(color = group_color[match(group, group_order)])
    p <- ggplot(plotdat_mean, aes(x = PC1, y = PC2, fill = group)) +
      geom_vline(xintercept = 0, color = "gray50", size = .3, lty = "dashed") + 
      geom_hline(yintercept = 0, color = "gray50", size = .3, lty = "dashed")
    for (i in 1:nrow(plotdat_mean)) {
      path <-  plotdat_mean[i, 2:5] %>% 
        unlist() %>% 
        as.numeric()
      var_color <-  plotdat_mean[i, 6] %>%
        unlist() %>% 
        as.character()
      p <- p + annotate(geom = "segment", x = path[1], y = path[2], 
                        xend = path[3], yend = path[4], color = var_color, lwd = .4)
    }
    p <- p + geom_point(aes(color = group), size = 1, show.legend = F) +
      geom_point(data = plotdat_mean %>% 
                   dplyr::select(group, X1mean, X2mean) %>% 
                   unique(), 
                 aes(x = X1mean, y = X2mean, color = group), 
                 size = 2, inherit.aes = F) +
      stat_ellipse(aes(color = group, fill = group), geom = 'polygon', 
                   level = .75, alpha = .05, lty = 2, lwd = .3, show.legend = F) +
      scale_color_manual(values = group_color) +
      scale_fill_manual(values = group_color) +
      labs(x = xlab, y = ylab,
           title = "Principal Components Analysis", 
           color = "Group")
  } else if (display_type == "dot") {
    p <- ggplot(plotdat, aes(x = PC1, y = PC2, color = group)) +
      geom_vline(xintercept = 0, lty = 2, lwd = .4) +
      geom_hline(yintercept = 0, lty = 2, lwd = .4) +
      geom_point(size = 1.5) +
      stat_ellipse(aes(fill = group), geom = 'polygon', level = .75, 
                   alpha = .05, lty = 2, lwd = .3, show.legend = F) +
      scale_color_manual(values = group_color) +
      scale_fill_manual(values = group_color) +
      labs(x = xlab, y = ylab,
           title = "Principal Components Analysis", 
           color = "Group")
  } else {stop("ERROR in show_type parameter .. dot or line.")}
  p <- p +
    theme_bw() +
    theme(axis.ticks = element_line(linewidth = .4, color = "black"),
          axis.title = element_text(size = 8, color = "black"),
          axis.text = element_text(size = 8, color = "black"),
          axis.line = element_blank(),
          plot.title = element_text(size = 10, color = "black"),
          plot.subtitle = element_text(size = 10, color = "black"),
          panel.background = element_blank(),
          panel.border = element_rect(linewidth = .4, color = "black"),
          panel.grid = element_blank(),
          legend.text = element_text(size = 8, color = "black"),
          legend.title = element_text(size = 8, color = "black"),
          aspect.ratio = 3/4)
  cat("  ggsave(file = \"PCA.pdf\", width = 6, height = 4.5)\n")
  return(p)
}