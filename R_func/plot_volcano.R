# Info ----
# encoding: utf-8
# author: Jinxin Meng
# e-mail: mengjx855@163.com
# created data：2024-03-28
# modified data: 2024-03-28
# version: 0.1

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

## plot_volcano ----
# diff, input a diff-like data.frame.
# group_pair 传入一个向量；例如c("up","none","down)
# 颜色是默认的，中间的是没有富集feature的颜色，如果改的话，只改开头和结尾的颜色就可
# if x scale is not symmetric, the xlim is used to limit the x scale, eg. c(5, 5)
# names parameter rename colnames of diff file, a vector contain pval, log2FC and enriched.
# xlab and ylab specify x and y axis title, respectively.
# pval and log2fc specify horizonal and vertical lines position.
plot_volcano <- function(dat, dat_colnames = NULL,
                         group_order = c("Up", "None", "Down"), group_color = NULL,
                         pval = 0.05, log2FC = 1, xlim = NULL, title = "Volcano plot",
                         xlab = expression(log[2]*FoldChange), ylab = expression(-log[10]*pval),
                         aspect.ratio = 3/4, subtitle_keywords = c("FoldChange", "pval")) {
  dat <- data.frame(dat, check.names = F)
  if (!all(c("name", "log2FC", "pval") %in% colnames(dat)) & is.null(dat_colnames)) stop("group field (name|log2FC|pval|..)")
  if (!is.null(dat_colnames)) dat <- data.frame(dat, check.names = F) %>% dplyr::rename(all_of(dat_colnames))
  if (! "enriched" %in% colnames(dat)) dat <- mutate(dat, enriched = ifelse(log2FC > log2FC & pval < pval, "Up", ifelse(log2FC < -log2FC & pval < pval, "Down", "None")))
  if (!is.null(xlim)) dat <- mutate(dat, log2FC = ifelse(log2FC > xlim[1], xlim[1], log2FC), log2FC = ifelse(log2FC < xlim[2], xlim[2], log2FC))
  if (is.null(group_order)) group_order <- unique(dat$enriched)
  if (is.null(group_color)) group_color <- structure(c("#f46d43", "grey75", "#66bd63"), names = group_order)
  if (!is.null(group_order) & !is.null(group_color)) group_color <- structure(group_color, names = group_order)
  
  dat <- mutate(dat, enriched = factor(enriched, group_order))
  
  group_label <- structure(names = group_order,
                           c(paste0(group_order[1], " (n=", sum(diff$enriched == group_order[1]), ")"),
                             paste0(group_order[2], " (n=", sum(diff$enriched == group_order[2]), ")"),
                             paste0(group_order[3], " (n=", sum(diff$enriched == group_order[3]), ")")))
  
  subtitle <- paste0("'Coefficient: ", subtitle_keywords[1], "'~'>'~'", 2^log2FC,
                     ",'~~italic('", subtitle_keywords[2], "')~'<'~'", pval, "'") %>% 
    as.formula() %>%
    eval()
  
  p <- ggplot(dat, aes(log2FC, -log10(pval), color = enriched)) +
    geom_point(size = .7) + 
    scale_color_manual(values = group_color, labels = group_label) + 
    scale_y_continuous(expand = c(.01, .01)) +
    labs(x = xlab, y = ylab, title = title, color = "Enriched in", subtitle = subtitle) + 
    geom_vline(xintercept = c(-log2FC, log2FC), lty = 2, lwd = .4, color = "black") + 
    geom_hline(yintercept = -log10(pval), lty = 2, lwd = .4, color = "black") +
    theme_bw() +
    theme(axis.ticks = element_line(linewidth = .4, color = "black"),
          axis.text = element_text(size = 8, color = "black"),
          axis.title = element_text(size = 8, color = "black"),
          axis.line = element_blank(),
          plot.title = element_text(size = 10, color = "black"),
          panel.background = element_rect(linewidth = .4, color = "black"),
          panel.grid = element_blank(),
          legend.position.inside = c(.02, .98),
          legend.justification = c(0, 1),
          legend.text = element_text(size = 8, color = "black", hjust = 1),
          legend.title = element_text(size = 8, color = "black"),
          legend.background = element_blank(),
          legend.key.height = unit(6, "mm"),
          legend.spacing.x = unit(.2, "mm"),
          aspect.ratio = aspect.ratio)
  
  message("  ggsave(file = \"diff_volcano_scatterplot.pdf\", width = 6, height = 4.5)\n")
  return(p)
}



