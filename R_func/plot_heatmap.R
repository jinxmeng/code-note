# Encoding: utf-8
# Author: Jinxin Meng
# Email: mengjx855@163.com
# Created Data：2024-03-09
# Modified Data: 2024-03-09
# Version: 0.1

library(tidyr)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
source("F:/Code/R_func/color.R")

#### plot_heatmap ####
# Jinxin Meng, 20240309
# profile format:
#       S1  S2  S3  S4
# OTU1  10  20  54  50
# OTU2  30  21  24  40
# OTU3  ..  ..  ..  ..
# row_annotation
# name  class ..
# OTU1  A     ..
# OTU2  B     ..
# ..    ..    ..
# col_annotation
# name  group
# S1    ctr
# S2    case
# ..    ..
# anno_colors
# list(group = c(ctr = "red", case = "black"),
#      class = c(A = "blue", B = "green"))
# show_row_anno, show_col_anno
# c("class", "group")

plot_heatmap <- function(profile, scale = "row", border_color = NA, title = "scale value",
                         cellwidth = 3, cellheight = 3,
                         cluster_rows = T, cluster_cols = T,
                         treeheight_row = 15, treeheight_col = 15,
                         show_rownames = F, show_colnames = F,
                         row_annotation = NULL, col_annotation  = NULL, annotation_legend = T, 
                         show_row_anno = NULL, show_col_anno = NULL, annotation_colors = NULL) {
  # 行列注释
  if (!is.null(row_annotation) & !(all(rownames(profile) %in% row_annotation$name))) stop("info. in row_annotation is lacking")
  if (!is.null(col_annotation) & !(all(colnames(profile) %in% col_annotation$name))) stop("info. in col_annotation is lacking")
  if (is.null(annotation_colors)) {
    annotation_colors <- list()
    if (!is.null(row_annotation) & !is.null(show_row_anno) & !all(show_row_anno %in% colnames(row_annotation))) stop("show_row_anno not in row_annotation")
    if (!is.null(row_annotation) & !is.null(show_row_anno)){
      row_annotation <- select(row_annotation, all_of(c("name", show_row_anno)))
      for (i in show_row_anno) {
        x = unique(row_annotation[,i])
        if (length(x) == 2) colors = structure(colors_pick("comparison2_2"), names = x)
        if (length(x) == 3) colors = structure(colors_pick("comparison_3"), names = x)
        if (length(x) >= 4 & length(x) <= 8) colors = structure(colors_pick("set2_8", length(x)), names = x)
        if (length(x) > 8 & length(x) <= 12) colors = structure(colors_pick("set3_12", length(x)), names = x)
        if (length(x) > 12 & length(x) <= 20) colors = structure(colors_pick("paired_20", length(x)), names = x)
        if (length(x) > 20) colors = structure(colors_continuous("hue", length(x)), names = x)
        annotation_colors[[i]] <- colors
      }
    }
    if (!is.null(col_annotation) & !is.null(show_col_anno) & !all(show_col_anno %in% colnames(col_annotation))) stop("show_col_anno not in col_annotation")
    if (!is.null(col_annotation) & !is.null(show_col_anno)) {
      col_annotation <- select(col_annotation, all_of(c("name", show_col_anno)))
      for (i in show_col_anno) {
        x = unique(col_annotation[,i])
        if (length(x) == 2) colors = structure(colors_pick("comparison2_2"), names = x)
        if (length(x) == 3) colors = structure(colors_pick("comparison_3"), names = x)
        if (length(x) >= 4 & length(x) <= 8) colors = structure(colors_pick("set2_8", length(x)), names = x)
        if (length(x) > 8 & length(x) <= 12) colors = structure(colors_pick("set3_12", length(x)), names = x)
        if (length(x) > 12 & length(x) <= 20) colors = structure(colors_pick("paired_20", length(x)), names = x)
        if (length(x) > 20) colors = structure(colors_continuous("hue", length(x)), names = x)
        annotation_colors[[i]] <- colors
      }
    }
  }
  row_annotation <- column_to_rownames(row_annotation, var = "name")
  col_annotation <- column_to_rownames(col_annotation, var = "name")
  pheatmap(profile, scale = scale,
          cluster_rows = cluster_rows, cluster_cols = cluster_cols,
          treeheight_row = treeheight_row, treeheight_col = treeheight_col,
          show_rownames = show_rownames, show_colnames = show_colnames,
          cellwidth = cellwidth, cellheight = cellwidth,
          border_gp = gpar(col = "black"), border_color = border_color,
          annotation_row = row_annotation, annotation_col = col_annotation,
          annotation_colors = annotation_colors,
          annotation_legend = annotation_legend,
          heatmap_legend_param = list(border = "black",  # 修改图例
                                      title = title, 
                                      title_gp = grid::gpar(fontface = "plain", fontsize = 10),
                                      title_position = "topleft",
                                      legend_direction = "vertical",
                                      legend_width = unit(4, "cm"),
                                      labels_gp = grid::gpar(fontsize = 8)))
}

