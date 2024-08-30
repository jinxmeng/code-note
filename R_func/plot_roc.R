# Encoding: utf-8
# Author: Jinxin Meng
# created date: 2023-8-16
# modified date: 2023-8-16

library(dplyr)
library(pROC)
library(ggplot2)

#### Function1 ####
# 针对一个roc对象绘制图，可输入曲线的颜色
plot_roc <- function(roc, color = "#238443", plot_se = F){
  label <- paste0("AUC: ", round(roc$auc, digits = 3), 
                  "\n(95% Cl: ", paste(round(ci.auc(roc), digits = 3)[c(1,3)], collapse = " ~ "), ")")
  p <- ggroc(roc, legacy.axes = T, lwd = .4, color = color) +
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color = "black", lty = "dashed", lwd = .2) +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    annotate("text" ,x = 0.75, y = 0.125 , label = label, size = 3) +
    # scale_x_continuous(expand = c(.005, .005)) +
    # scale_y_continuous(expand = c(.005, .005)) +
    theme_bw() +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.ticks = element_line(linewidth = .5, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.line = element_blank(),
          panel.background = element_rect(linewidth = .4, color = "black"),
          panel.grid = element_blank(),
          aspect.ratio = 1)
  if(isTRUE(plot_se)) {
    roc_se <- ci.se(roc, specificities = seq(0, 1, 0.01), conf.level = 0.95) %>% 
      data.frame(check.names = F) %>% 
      dplyr::rename(lower = all_of("2.5%"), upper = all_of("97.5%"), median = all_of("50%")) %>% 
      rownames_to_column(var = "spec") %>% 
      mutate(spec = as.numeric(spec))
    p <- p + geom_ribbon(data = roc_se, aes(x = 1 - spec, ymin = lower, ymax = upper), fill = color, alpha = .1)
  }
  cat("  ggsave(file = \"roc.pdf\", width = 4, height = 4)\n")
  return(p)
} 

#### Function2 ####
# 针对一个roc列表对象绘制图，可输入曲线的颜色，用向量
plot_roc_list <- function(roc_list, colors = NULL, plot_se = F){
  if (is.null(colors)) {
    colors <- c("#1f78b4","#33a02c","#e31a1c","#ff7f00","#6a3d9a", 
                "#a6cee3","#b2df8a","#fb9a99","#fdbf6f","#cab2d6")
    colors <- colors[1:length(roc_list)]
  } else { colors <- colors }
  labels <- lapply(roc_list, \(x) 
                   paste0(round(x$auc, digits = 3), "\n(95% Cl: ", paste(round(ci.auc(x), digits = 3)[c(1,3)], collapse = " ~ "), ")")) %>% 
    unlist() %>% paste0(names(roc_list), " AUC: ", .)
  p <- ggroc(roc_list, legacy.axes = T, lwd = .4) +
      geom_segment(data = data.frame(x = 0, y = 0), aes(x = x, y = y, xend = 1, yend = 1),
                   color = "black", lty = "dashed", lwd = .4) +
      scale_color_manual(values = colors, breaks = names(roc_list), labels = labels) +
      scale_x_continuous(expand = c(.005, .005)) +
      scale_y_continuous(expand = c(.005, .005)) +
      labs(x = "1 - Specificity", y = "Sensitivity", color = "") +
      theme_bw() +
      theme(axis.text = element_text(size = 8, color = "black"),
            axis.ticks = element_line(linewidth = .4, color = "black"),
            axis.title = element_text(size = 10, color = "black"),
            axis.line = element_blank(),
            panel.background = element_rect(linewidth = .4, color = "black"),
            panel.grid = element_blank(),
            legend.position = c(.95, .02),
            legend.justification = c(1, 0),
            legend.text = element_text(size = 8, color = "black"),
            legend.background = element_blank(),
            legend.key.height = unit(7, "mm"),
            legend.spacing.x = unit(.8, "mm"),
            aspect.ratio = 1)
  if (isTRUE(plot_se)) {
    roc_se <- lapply(roc_list, \(x) ci.se(x, specificities = seq(0, 1, 0.01), conf.level = 0.95) %>% 
                       data.frame(check.names = F) %>% 
                       rename(lower = all_of("2.5%"), upper = all_of("97.5%"), median = all_of("50%")) %>% 
                       rownames_to_column(var = "spec") %>% 
                       mutate(spec = as.numeric(spec))) %>% 
      purrr::map2_dfr(., names(.), \(x, y) add_column(x, class = y))
    p <- p + geom_ribbon(data = roc_se, aes(x = 1 - spec, ymin = lower, ymax = upper, fill = class), alpha = .1, 
                         inherit.aes = F, show.legend = F) +
      scale_fill_manual(values = colors)
  }
  cat("  ggsave(file = \"roc_list.pdf\", width = 4, height = 4)\n")
  return(p)
}
