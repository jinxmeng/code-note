# Encoding: utf-8
# Author: Jinxin Meng
# Created Data：2022-9-27
# Modified Date: 2022-9-28
# Version: 0.5

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(vegan)

#### Function1 dbRDA ####
# otu输入正常格式的otu
# group输入一个group表，目前支持一个类型的分组，还没有加入环境矩阵
# group_order输入分组的顺序
# group_color输入分组的颜色
# 此函数计算了bary-curtis距离的dbRDA，输出的RDA为I型标尺
# 此函数计算了dbRDA的R2adj, 并进行PERMANOVA检验，在图的subtitle展示
# 输出为ggplot对象
plot_dbRDA <- function(otu, group, group_order = NULL, group_color = NULL) {
  # group的level
  if(is.null(group_order)){
    group_order <- unique(group$group)
  } else {
    group_order <- group_order
  }
  
  # group的配色
  if(is.null(group_color)){
    color <- c("#1f78b4","#33a02c","#e31a1c","#ff7f00","#6a3d9a","#ffff99","#b15928","#a6cee3","#b2df8a","#fb9a99","#fdbf6f","#cab2d6")
    group_color <- rep(color, times = ceiling(length(group_order)/length(color)))[seq_len(length(group_order))]
  } else {
    group_color <- group_color
  }
  
  group2 <- group[match(colnames(otu), group$sample),] %>% 
    data.frame(check.names = F, row.names = NULL) %>% 
    column_to_rownames(var = "sample")
  
  # 基于距离的矩阵
  distance <- vegdist(t(otu), method = "bray")
  
  # db-RDA
  dbRDA <- capscale(distance ~ ., group2)
  # dbRDA2 <- dbrda(distance ~ ., group)
  dbRDA_summary <- summary(dbRDA, scaling = 1) # I型标尺
  
  # 样方坐标
  dbRDA_point <- dbRDA_summary$sites %>%
    data.frame(check.names = F) %>% 
    select(all_of(colnames(.)[1:2])) %>% 
    rownames_to_column(var = "sample") %>% 
    mutate(group = group$group[match(sample, group$sample)])
  
  # 效应变量坐标
  dbRDA_variable <- dbRDA_summary$centroids %>% 
    data.frame(check.names = F) %>% 
    select(all_of(colnames(.)[1:2])) %>% 
    rownames_to_column(var = "group") %>% 
    mutate(group = gsub("group", "", group))
  
  # 前两轴的解释度
  dbRDA_eig <- round(dbRDA_summary$cont$importance[2,]*100, digits = 2)
  
  # R2矫正
  R2adj <- round(RsquareAdj(dbRDA)$adj.r.squared*100, digits = 4)
  label <- paste0("db-RDA Constrained: R2adj = ", R2adj, "%")

  # PERMANOVA
  group <- group[match(rownames(as.matrix(distance)), group$sample),] %>% 
    as.data.frame(row.names = NULL)  # 整体水平比较,保证样本在两个数据集中对应
  adonis <- adonis2(distance ~ group, group, permutations = 999)  # 基于距离矩阵计算
  adonis_p <- adonis[1,5]  # 提取p
  adonis_R2 <- round(adonis[1,3]*100, digits = 4)  # 提取R2
  label2 <- paste0("PERMANOVA: R2 = ",adonis_R2,"%  p = ", adonis_p)
  
  x <- colnames(dbRDA_point)[2]
  y <- colnames(dbRDA_point)[3]
  
  dbRDA_point$group <- factor(dbRDA_point$group, levels = group_order)
  dbRDA_point$group <- factor(dbRDA_point$group, levels = group_order)
  
  p <- ggplot() +
    geom_vline(xintercept = 0, color = "gray50", size = .3, lty = "dashed") + 
    geom_hline(yintercept = 0, color = "gray50", size = .3, lty = "dashed") +
    geom_point(data = dbRDA_point, aes_string(x = x, y = y, color = "group"), size = 1) +
    stat_ellipse(data = dbRDA_point, aes_string(x = x, y = y, fill = "group", color = "group"), 
                 geom = 'polygon', alpha = .02, level = .95, show.legend = F) +
    geom_segment(data = dbRDA_variable, aes_string(x = 0, y = 0, xend = x, yend = y), 
                 arrow = arrow(length = unit(.8, 'mm')), size = .3, color = "black") +
    ggrepel::geom_text_repel(data = dbRDA_variable, aes_string(x = x , y = y, label = "group"), 
                             color = 'black', size = 2) +
    scale_color_manual(values = group_color) +
    scale_fill_manual(values = group_color) +
    labs(x = paste0(x, " (", dbRDA_eig[1], "%)"), y = paste0(y, " (", dbRDA_eig[2], "%)"),
         title = "bray-curtis dbRDA", subtitle = paste0(label,"\n", label2), color = "Group") +
    theme_bw() +
    theme(axis.ticks = element_line(linewidth = .4, color = "black"),
          axis.title = element_text(size = 8, color = "black"),
          axis.text = element_text(size = 8, color = "black"),
          axis.line = element_blank(),
          plot.title = element_text(size = 10, color = "black"),
          plot.subtitle = element_text(size = 8, color = "black"),
          panel.background = element_blank(),
          panel.border = element_rect(linewidth = .4, color = "black"),
          panel.grid = element_blank(),
          legend.text = element_text(size = 8, color = "black"),
          legend.title = element_text(size = 8, color = "black"),
          aspect.ratio = 3/4)
  cat("Suggestion: ggsave(\"dbRDA.pdf\", width = 6, height = 4.5)")
  return(p)
}