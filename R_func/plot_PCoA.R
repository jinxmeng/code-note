#### info ####
# encoding: utf-8
# author: Jinxin Meng
# e-mail: mengjx855@163.com
# created data：2022-05-29
# modified data: 2024-03-04
# version: 0.5

# 2022-10-01: 添加选择不同距离尺度的参数 dis_method
# 2023-01-01: update function
# 2023-12-04: update function: check_file_name was deprecated.

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(vegan)

#### calcu_PCoA ####
# PCoA analysis
# input: profile | distance
# dim: how many dimension is used to analyze
# cumulative_eig: 轴累积解释度达到多少阈值，输出多少个轴
# prefix: 修改列名的前缀
# add_eig: 列名中含有解释度
# dis_method: 仅支持vegan包中的距离计算方法
calcu_PCoA <- function(profile, distance, group, group_colnames = NULL, dim = 2, cumulative_eig = NULL, dis_method = "bray", 
                       prefix = NULL, adonis2 = T) {
  if (!missing(profile) & missing(distance)) {
    distance <- vegdist(t(data.frame(profile, check.names = F)), method = dis_method, na.rm = T)
  } else if (missing(profile) & !missing(distance)) {
    distance <- distance
  } else {
    stop("Error in distance compute: cannot open otu or distance data.")
  }
  
  if (is.null(cumulative_eig)) {
    PCoA <- cmdscale(distance, k = dim, eig = T)
    PCoA_eig <- round(PCoA$eig/sum(PCoA$eig) * 100, digits = 2)
    PCoA_points <- rownames_to_column(data.frame(PCoA$points), var = "sample")
  } else if (!is.null(cumulative_eig)) {
    dim = ncol(profile) - 1
    PCoA <- cmdscale(distance, k = dim, eig = T) %>% suppressWarnings()
    PCoA_eig <- round(PCoA$eig/sum(PCoA$eig) * 100, digits = 2)
    tmp_vec <- cumsum(PCoA_eig)
    for (i in 1:length(tmp_vec)) {
      if (tmp_vec[i] >= cumulative_eig) {
        dim <- i
        break
      }
    }
    PCoA_points <-  dplyr::select(data.frame(PCoA$points), all_of(1:dim)) %>% 
      rownames_to_column(var = "sample")
  }
  
  if (!is.null(prefix)) {
    prefix <- as.character(prefix)
    colnames(PCoA_points)[2:(dim + 1)] <- paste0(prefix, "_", seq(1, dim))
  } else {
    colnames(PCoA_points)[2:(dim + 1)] <- paste0("PCoA", seq(1, dim))
  }
  
  out = list(dist = distance, 
             points = PCoA_points, 
             group = group,
             dim = dim, 
             eig = PCoA_eig, 
             eig_ = paste0(colnames(PCoA_points)[2:(dim + 1)], "(", PCoA_eig[1:dim], "%)"))
  
  if (isTRUE(adonis2)) {
    if (missing(group)) stop("need group file when perform adonis2 analsis")
    if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_colnames)) stop("group field (sample|group)")
    if (!is.null(group_colnames)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_colnames))
    group <- group[match(rownames(as.matrix(distance)), group$sample),] %>% 
      as.data.frame(row.names = NULL)
    adonis <- adonis2(distance ~ group, group, permutations = 999)
    label <- paste0("R2=", round(adonis[1,3], digits = 4), " p<", adonis[1,5])
    
    calcu_adjusted_r2 <- function(adonis_object) {
      n_observations <- adonis_object$Df[3]+1
      d_freedom <- adonis_object$Df[1]
      r2 <- adonis_object$R2[1]
      adjusted_r2 <- RsquareAdj(r2, n_observations, d_freedom)
      return(adjusted_r2)
    } 
    
    r2adj <- calcu_adjusted_r2(adonis)
    out[["group"]] = group
    out[["adonis2"]] = adonis
    out[["adonis2_r2"]] = round(adonis[1,3], digits = 4)
    out[["adonis2_r2adj"]] = round(r2adj, digits = 4)
    out[["adonis2_p"]] = adonis[1,5]
    out[["label"]] = label
  }
  return(out)
}

#### plot_PCoA ####
# PCoA analysis
# input: profile & sample_group
# group_colnames: a vector for current colnames of the data.frame with new name.
# for example: group_colnames = c(sample = "Sample", val = "Val")
# group_order: group order
# group_color: group color
# dis_method: 仅支持vegan包中的距离计算方法
plot_PCoA <- function(profile, group, distance, 
                      group_order = NULL, group_color = NULL, group_colnames = NULL,
                      dis_method = "bray", display_type = "line", ellipse_level = .75, 
                      title = NULL, add_lab_to_plot = F, show_legend = T) {
  if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_colnames)) stop("group field (sample|group)")
  if (!is.null(group_colnames)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_colnames))
  if (is.null(group_order)) group_order <- unique(group$group)
  if (is.null(group_color)) {
    color <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494","#b3b3b3")
    group_color <- rep(color, times = ceiling(length(group_order)/length(color)))[1:length(group_order)]
  }
  if (is.null(title)) title = paste0(stringr::str_to_sentence(string = dis_method), "-distance PCoA")
  
  if (!missing(profile) & missing(distance)) {
    distance <- vegdist(t(data.frame(profile, check.names = F)), method = dis_method, na.rm = T)
  } else if (missing(profile) & !missing(distance)) {
    distance <- distance
  } else {
    stop("Error in distance compute: cannot open otu or distance data.")
  }
  
  PCoA <- cmdscale(distance, k = 2, eig = T)
  PCoA_points <- rownames_to_column(data.frame(PCoA$points), var = "sample")
  PCoA_eig <- round(PCoA$eig/sum(PCoA$eig) * 100, digits = 2)
  plotdat <- merge(PCoA_points, group %>% select(sample, group), by = "sample", all.x = T) %>% 
    mutate(group = factor(group, group_order))
  
  # adonis
  group <- group[match(rownames(as.matrix(distance)), group$sample),] %>%
    as.data.frame(row.names = NULL)
  adonis <- adonis2(distance ~ group, group, permutations = 999)
  # label <- paste0("R2 = ",round(adonis[1,3], digits = 4),"  p = ", adonis[1,5])
  # label <-'R'^2~'='~'0.4226'~~italic('p')~'<'~'0.001'
  label <- paste0("'R'^2~'='~'", round(adonis[1,3], digits = 4), "'~~italic('p')~'<'~'", adonis[1,5], "'") %>% 
    as.formula() %>% eval()
  
  if (display_type == "line") {
    plotdat_mean <- plotdat %>% select(-sample) %>% group_by(group) %>% 
      summarise_all(mean) %>% rename(X1mean = X1, X2mean = X2) %>% 
      merge(x = plotdat %>% select(-sample) %>% relocate(group), y = ., by = "group") %>% 
      mutate(color = group_color[match(group, group_order)]) 
    
    p <- ggplot(plotdat_mean, aes(x = X1, y = X2, fill = group)) +
      geom_vline(xintercept = 0, lty = 2, lwd = .4) +
      geom_hline(yintercept = 0, lty = 2, lwd = .4)
    
    for (i in 1:nrow(plotdat_mean)) {
      path <- plotdat_mean[i, 2:5] %>% unlist() %>% as.numeric()
      var_color <-  plotdat_mean[i, 6] %>% unlist() %>% as.character()
      p <- p + annotate(geom = "segment", x = path[1], y = path[2], 
                        xend = path[3], yend = path[4], color = var_color, lwd = .4)
    }
    
    p <- p + geom_point(aes(color = group), size = 1, show.legend = F) +
      geom_point(data = plotdat_mean %>% select(group, X1mean, X2mean) %>% unique(), 
                 aes(x = X1mean, y = X2mean, color = group), size = 2, inherit.aes = F) +
      stat_ellipse(aes(color = group, fill = group), geom = 'polygon', level = ellipse_level, 
                   alpha = .05, lty = 2, lwd = .3, show.legend = F) +
      scale_color_manual(values = group_color) +
      scale_fill_manual(values = group_color) +
      labs(x = paste("PCoA1 (", PCoA_eig[1], "%)", sep = ""), 
           y = paste("PCoA2 (", PCoA_eig[2], "%)", sep = ""),
           title = title, color = "Group", subtitle = label) +
      theme_bw() +
      theme(axis.ticks = element_line(linewidth = .4, color = "black"),
            axis.title = element_text(size = 8, color = "black"),
            axis.text = element_text(size = 8, color = "black"),
            axis.line = element_blank(),
            plot.title = element_text(size = 10, color = "black"),
            plot.subtitle = element_text(size = 10, color = "black"),
            panel.border = element_rect(linewidth = .4, color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_line(linewidth = .4, color = "grey90"),
            panel.grid.minor = element_line(linewidth = .4, color = "grey90"),
            legend.text = element_text(size = 8, color = "black"),
            legend.title = element_text(size = 8, color = "black"),
            aspect.ratio = 3/4)
    
    if (isTRUE(add_lab_to_plot)) p <- p + ggrepel::geom_label_repel(
      data = unique(select(plotdat_mean, group, X1mean, X2mean)),
      aes(x = X1mean, y = X2mean, label = group), size = 1.5, show.legend = F)
    
    if (isFALSE(show_legend)) p = p + guides(color = "none")
  
  } else if (display_type == "dot") {
    p <- ggplot(plotdat, aes(x = X1, y = X2, color = group)) +
      geom_vline(xintercept = 0, lty = 2, lwd = .4) +
      geom_hline(yintercept = 0, lty = 2, lwd = .4) +
      geom_point(size = 1) +
      stat_ellipse(aes(fill = group), geom = 'polygon', level = .75, 
                   alpha = .05, lty = 2, lwd = .3, show.legend = F) +
      scale_color_manual(values = group_color) +
      scale_fill_manual(values = group_color) +
      labs(x = paste("PCoA1 (", PCoA_eig[1], "%)", sep = ""), 
           y = paste("PCoA2 (", PCoA_eig[2], "%)", sep = ""),
           title = title, color = "Group", subtitle = label) +
      theme_bw() +
      theme(axis.ticks = element_line(linewidth = .4, color = "black"),
            axis.title = element_text(size = 8, color = "black"),
            axis.text = element_text(size = 8, color = "black"),
            axis.line = element_blank(),
            plot.title = element_text(size = 10, color = "black"),
            plot.subtitle = element_text(size = 10, color = "black"),
            panel.border = element_rect(linewidth = .4, color = "black"),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            legend.text = element_text(size = 8, color = "black"),
            legend.title = element_text(size = 8, color = "black"),
            aspect.ratio = 3/4)
  } else {
    stop("ERROR in show_type parameter .. dot or line.")
  }
  message("  ggsave(file = \"div_PCoA.pdf\", width = 6, height = 4.5)")
  return(p)
}

#### calcu_pairwise_adonis ####
# PREMANOVA analysis
# input: profile & mapping
# group_order: group order
calcu_pairwise_adonis <- function(profile, group, 
                                  group_order = NULL,group_names = NULL, group_colnames = NULL,
                                 dis_method = "bray", permutations = 999) {
  profile <- data.frame(profile, check.names = F)
  if (!is.null(group_colnames)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_colnames))
  if (is.null(group_order)) group_order <- unique(group$group)
  
  adonis2 <- data.frame(group_pair = character(), r2 = numeric(), r2adj = numeric(), pval = numeric())
  
  calcu_adjusted_r2 <- function(adonis_object) {
    n_observations <- adonis_object$Df[3]+1
    d_freedom <- adonis_object$Df[1]
    r2 <- adonis_object$R2[1]
    adjusted_r2 <- RsquareAdj(r2, n_observations, d_freedom)
    return(adjusted_r2)
  }
  
  for (i in 1:(length(group_order) - 1)) {
    for (j in (i + 1):length(group_order)) {
      group_ij <- subset(group, group %in% c(group_order[i], group_order[j]))
      otu_ij <- profile[, group_ij$sample]
      adonis_ij <- adonis2(t(otu_ij) ~ group, group_ij, permutations = permutations, distance = dis_method)
      r2adj <- calcu_adjusted_r2(adonis_ij)
      adonis2 <- adonis2 %>% add_row(group_pair = paste0(group_order[i],'_vs_',group_order[j]),
                                     r2 = adonis_ij$R2[1], r2adj = r2adj, pval = adonis_ij$`Pr(>F)`[1])
    }
  }
  return(adonis2)
}


#### calcu_distance ####
calcu_distance <- function(profile, method = "unifrac", tree, weighted = T, ...) {
  if (method == "unifrac") {
    if (missing(tree)) stop("need phylogentic tree if calcu unifrac-based distance\n🥰 input a tree: ape::read.tree(\"genomospecies.tre\")")
    ps <- phyloseq::phyloseq(phyloseq::otu_table(profile, taxa_are_rows = T), tree)
    distance <- phyloseq::UniFrac(ps, weighted = weighted)
    if (isTRUE(weighted)) message("Output: weighted-unifrac-based distance.")
    if (isFALSE(weighted)) message("Output: unweighted-unifrac-based distance.")
  } else {
    distance <- vegdist(t(data.frame(profile, check.names = F)), method = method, na.rm = T)
    message(paste0("Output: ", method ,"-based distance."))
  }
  return(distance)
}


