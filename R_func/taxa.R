#### info ###
# encoding: utf-8
# author: Jinxin Meng
# e-mail: mengjx855@163.com
# created data：2022-9-18
# modified date: 2024-03-05
# version: 0.1

library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
source("F:/code/r_func/profile_process.R")

# taxa_split ####
# 拆分物种信息表
# taxonomy: a data.frame containning two field, name|taxonomy..., also rename by taxonomy_name.
# sep: the separate label.
# na_fill: unclassified rename.
# rm_suffix: remove taxa name with suffix providing by GTDB.
# taxa_level = "d:s", 指定层级列表
taxa_split <- function(taxonomy, taxonomy_colnames = NULL, sep = ";", taxa_level = "d:s", 
                       na_fill = "Unknown", rm_suffix = F) {
  if (!all(c("name", "taxonomy") %in% colnames(taxonomy)) & is.null(taxonomy_colnames)) stop("taxonomy field (name|taxonomy)")
  if (!is.null(taxonomy_colnames)) taxonomy <- dplyr::rename(data.frame(taxonomy, check.names = F), all_of(taxonomy_colnames))

  dat <- data.frame(taxonomy, check.names = F)
  # 判断层级
  taxa_level_vec <- c("d", "p", "c", "o", "f", "g", "s", "t")
  taxa_name_vec <- c("domain", "phylum", "class", "order", "family", "genus", "species", "strain")
  cut_off <- c(match(unlist(strsplit(taxa_level, ":"))[1], taxa_level_vec):match(unlist(strsplit(taxa_level, ":"))[2], taxa_level_vec))
  taxa_name_vec <- taxa_name_vec[cut_off]
  taxa_level_vec <- taxa_level_vec[cut_off] %>% paste0(., "__")
  # 处理
  res <- map_dfr(dat$taxonomy, \(x) 
      if(is.na(x)){ purrr::map_vec(taxa_level_vec, \(i) ifelse(grepl("__$", i), paste0(i, na_fill), i)) %>% t() %>% data.frame() } 
      else { purrr::map_vec(unlist(strsplit(x, sep)), \(i) ifelse(grepl("__$", i), paste0(i, na_fill), i))[cut_off] %>% t() %>% data.frame() })
  colnames(res) <- taxa_name_vec
  res <- add_column(res, name = dat$name, .before = 1)
  # 去除门的子群标记，例如，p__Firmicutes_A, 修改为p__Firmicutes
  if (isTRUE(rm_suffix)) {
    res <- apply(res, 2, \(x) gsub("_\\w$", "", x = x) %>% gsub("_\\w ", " ", x = ., perl = T)) %>% data.frame(check.names = F)
  }
  return(res)
}

# taxa_trans ####
# 转换物种表的层级关系
# otu: input a otu table
# taxonomy: input a data.frame, colnames following: name|phylum|family..., also rename by taxonomy_name parameter.
# taxonomy_name = a vector used to rename taxon table.
# to: 指定转变为哪个级别的
# top_n: select top n name, including Other.
# top_list: select top n name by specify name.
# other_name: a name to rename other name.
# out_all: output all result.
# na_fill: if some name a not corresponding front level, rename it.
# smp2grp: merge sample to group.
# group: input a data.frame, colnames following: sample|group, also rename by group_name parameter.
# method: merge method, normal by mean or median.
# reture a otu table.
taxa_trans <- function(otu, taxonomy, group, to = "family", top_n = 12, top_list = NULL, 
                       other_name = "Other", out_all = F, na_fill = "Unclassified",
                       transRA = F, smp2grp = F, method = "mean", 
                       taxonomy_colnames = NULL, group_colnames = NULL) {
  otu <- data.frame(otu, check.names = F)
  taxonomy <- data.frame(taxonomy, check.names = F)
  if (!all(c("name", to) %in% colnames(taxonomy)) & is.null(taxonomy_colnames)) stop(paste0("taxonomy field must have \"name\" and \"", to, "\", and other level optional input."))
  if (!is.null(taxonomy_colnames)) taxonomy <- dplyr::rename(taxonomy, all_of(taxonomy_colnames))
  if(isTRUE(out_all)) top_n <- 0
  dat <- otu %>% 
    mutate(taxa = taxonomy[match(rownames(.), taxonomy$name), to],
           taxa = ifelse(is.na(taxa), na_fill, taxa)) %>% 
    group_by(taxa) %>% summarise_all(sum) %>% ungroup() %>% 
    column_to_rownames(var = "taxa") %>% data.frame(check.names = F)
  # trans taxon
  if(top_n > 0 & is.null(top_list)) {
    taxa_vec <- data.frame(val = rowSums(dat)) %>%
      arrange(desc(val)) %>% head(n = top_n - 1) %>% 
      rownames(.) %>% unlist() %>% as.character()
    dat <- dat %>% rownames_to_column(var = "name") %>% 
      mutate(name = ifelse(name %in% taxa_vec, name, other_name)) %>% 
      group_by(name) %>% summarise_all(sum) %>% ungroup() %>% 
      data.frame(check.names = F) %>% column_to_rownames(var = "name")
    } else if(!is.null(top_list)) {
      dat <- dat %>% rownames_to_column(var = "name") %>% 
        mutate(name = ifelse(name %in% top_list, name, other_name)) %>% 
        group_by(name) %>% summarise_all(sum) %>% ungroup() %>% 
        data.frame(check.names = F) %>% column_to_rownames(var = "name")
    }
  # sample to group
  if (isTRUE(smp2grp)) {
    if (missing(group)) stop("missing group file.")
    if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_colnames)) stop("group field (sample|group)")
    if (!is.null(group_colnames)) group <- dplyr::rename(data.frame(group, check.names = F), all_of(group_colnames))
    dat <- data.frame(t(dat), check.names = F) %>% 
      mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
      group_by(group) %>% summarise_all(method) %>% ungroup() %>% 
      column_to_rownames(var = "group") %>% t() %>% data.frame(check.names = F)
  }
  # transRA
  if(isTRUE(transRA)) {
    dat <- apply(dat, 2, \(x) x/sum(x)*100) %>% data.frame(check.names = F)
  }
  dat <- dat[rowSums(dat)!=0, colSums(dat)!=0]
  return(dat)
}

# plot_compos ####
# 物种组成,
# otu表输入行名为name，列名为样本的表
# group输入列名为sample和group的表
# taxonomy输入taxa.split()函数输出的表格，或者是指定新表格，第一列为name，随后是不同层级的分类信息
# display指定样本级别的组成还是分组级别的组成， 必须指定为sample或者group
# taxa_level显示哪个水平组成, 可选phylum, class, family, genus, species
# top_n输入显示多少个物种，不足数量的按最大数量算
# top_list指定物种列表
# group_order/sample_order，分别之指定sample和group的顺序
# taxa_order指定物种的排序
# taxa_color指定物种的颜色
# title指定图的标题内容
plot_compos <- function(otu, taxonomy, group, display = "group", taxa_level = "family", 
                        top_n = 12, top_list = NULL, group_order = NULL, sample_order = NULL,
                        taxa_order = NULL, taxa_color = NULL, plot_title = NULL,
                        taxonomy_colnames = NULL, group_colnames = NULL){
  if(missing(otu) & missing(taxonomy)) stop("missing parameter.")
  if(display == "group" & !missing(group)) {
    group <- data.frame(group, check.names = F) 
  } else if (display == "group" & missing(group)) { 
    stop("missing parameter.") 
  }
  # 读入数据
  otu <- data.frame(otu, check.names = F)
  taxonomy <- data.frame(taxonomy, check.names = F)
  colors <- c("#4E79A7FF","#A0CBE8FF","#F28E2BFF","#FFBE7DFF","#59A14FFF",
              "#8CD17DFF","#B6992DFF","#F1CE63FF","#499894FF","#86BCB6FF",
              "#E15759FF","#FF9D9AFF","#79706EFF","#BAB0ACFF","#D37295FF",
              "#FABFD2FF","#B07AA1FF","#D4A6C8FF","#9D7660FF","#D7B5A6FF")
  if(display == "group") {
    other_name <- paste0(str_to_lower(str_sub(taxa_level, 1, 1)), "__Other")
    dat <- taxa_trans(otu, taxonomy, group, to = taxa_level, top_n = top_n,
                      top_list = top_list, other_name = other_name, smp2grp = T, transRA = T,
                      group_colnames = group_colnames, taxonomy_colnames = taxonomy_colnames)
    if(is.null(taxa_order)) taxa_order <- names(sort(rowSums(dat), decreasing = T))
    if(is.null(taxa_color)) taxa_color <- rep(colors, time = ceiling(nrow(dat)/20))[1:nrow(dat)]
    if(is.null(group_order)) group_order <- names(sort(unlist(dat[taxa_order[1],])))
    if(is.null(plot_title)) plot_title = stringr::str_to_sentence(taxa_level)
    fill_title = paste0(plot_title, " taxa")
    plotdat <- rownames_to_column(dat, var = "name") %>%  
      gather(., key = "group", value = "val", -name) %>% 
      mutate(group = factor(group, group_order), 
             name = factor(name, taxa_order))
    p <- ggplot(plotdat, aes(group, val, fill = name)) +
      geom_bar(stat = "identity", position = position_stack(), color = "#000000", linewidth = .2, width = .8) +
      scale_fill_manual(values = taxa_color) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "", y = "Relative Abundance (%)", title = plot_title, fill = fill_title) +
      theme_classic() +
      theme(axis.line = element_line(linewidth = .4, color = "#000000"),
            axis.ticks = element_line(linewidth = .4, color = "#000000"),
            axis.text = element_text(size = 8, color = "#000000"),
            axis.title = element_text(size = 8, color = "#000000"),
            plot.title = element_text(size = 10, color = "#000000"),
            legend.text = element_text(size = 8, color = "#000000", face = "italic"),
            legend.title = element_text(size = 10, color = "#000000"),
            panel.grid = element_blank())
    message("  ggsave(file = \"compos_stacked_barplot.pdf\", width = 6, height = 4)")
  } else if(display == "sample") {
    other_name <- paste0(str_to_lower(str_sub(taxa_level, 1, 1)), "__Other")
    dat <- taxa_trans(otu, taxonomy, group, to = taxa_level, top_n = top_n,
                      top_list = top_list, other_name = other_name, smp2grp = F, transRA = T,
                      group_colnames = group_colnames, taxonomy_colnames = taxonomy_colnames)
    if(is.null(taxa_order)) taxa_order <- names(sort(rowSums(dat), decreasing = T))
    if(is.null(taxa_color)) taxa_color <- rep(colors, time = ceiling(nrow(dat)/20))[1:nrow(dat)]
    if(is.null(sample_order)) sample_order <- names(sort(unlist(dat[taxa_order[1],])))
    if(is.null(plot_title)) plot_title = stringr::str_to_sentence(taxa_level)
    fill_title = paste0(plot_title, " taxa")
    plotdat <- rownames_to_column(dat, var = "name") %>%  
      gather(., key = "sample", value = "val", -name) %>% 
      mutate(sample = factor(sample, sample_order),
             name = factor(name, taxa_order))
    p <- ggplot(plotdat, aes(sample, val, fill = name)) +
      geom_bar(stat = "identity", position = position_stack(), color = "#000000", linewidth = .2, width = 1) +
      scale_fill_manual(values = taxa_color) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "", y = "Relative Abundance (%)", title = plot_title, fill = fill_title) +
      theme_classic() +
      theme(axis.line = element_line(linewidth = .4, color = "#000000"),
            axis.ticks = element_line(linewidth = .4, color = "#000000"),
            axis.ticks.x = element_blank(),
            axis.text = element_text(size = 8, color = "#000000"),
            axis.text.x = element_blank(),
            axis.title = element_text(size = 8, color = "#000000"),
            plot.title = element_text(size = 10, color = "#000000"),
            legend.text = element_text(size = 8, color = "#000000", face = "italic"),
            legend.title = element_text(size = 10, color = "#000000"),
            panel.grid = element_blank())
    message("  ggsave(file = \"compos_stacked_barplot.pdf\", width = 8, height = 4)")
  } else {
    stop("error parameter.")
  }
  return(p)
}

# plot_mcompos ####
# 物种组成,
# otu表输入行名为name，列名为样本的表
# group输入列名为sample和group的表
# taxonomy输入taxa.split()函数输出的表格，或者是指定新表格，第一列为name，随后是不同层级的分类信息
# display指定样本级别的组成还是分组级别的组成， 必须指定为sample或者group
# top_n输入显示多少个分类，不足数量的按最大数量算
plot_mcompos <- function(otu, taxonomy, group, display = "group", top_n = 12, 
                         top_list = NULL, group_order = NULL, sample_order = NULL,
                         taxonomy_colnames = NULL, group_colnames = NULL) {
  if(missing(otu) & missing(taxonomy)) stop("missing parameter.")
  taxa_level <- setdiff(colnames(taxonomy), c("name", "domain"))
  
  p_list <- list()
  for (i in taxa_level) {
    p <- plot_compos(otu = otu, taxonomy = taxonomy, group = group, display = display,
                     top_n = top_n, top_list = top_list, taxa_level = i,
                     group_order = group_order,sample_order = sample_order, 
                     taxonomy_colnames = taxonomy_colnames, group_colnames = group_colnames) %>% 
      suppressMessages()
    p_list[[i]] <- p
  }
  p <- cowplot::plot_grid(plotlist = p_list, nrow = 2, align = "v")
  if (display == "group") {
    cat("  ggsave(file = \"compos_stacked_barplot.pdf\", width = 18, height = 8)")
  } else if (display == "sample") {
    cat("  ggsave(file = \"compos_stacked_barplot.pdf\", width = 30, height = 8)")
  }
  return(p)
}

# plot_compos_manual ----

plot_compos_manual <- function(otu, group, display = "group", top_n = 12, group_order = NULL, 
                               sample_order = NULL, taxa_order = NULL, taxa_color = NULL, plot_title = NULL,
                               group_colnames = NULL){
  if(display == "group" & !missing(group)) {
    group <- data.frame(group, check.names = F) 
  } else if (display == "group" & missing(group)) { 
    stop("missing parameter.") 
  }
  
  # otu <- apply(otu, 2, \(x) ifelse(is.nan(x), 0, x)) %>% 
  #   apply(., 2, \(x) ifelse(is.na(x), 0, x))

  otu <- otu[rowSums(otu)!=0, colSums(otu)!=0] %>% 
    data.frame(check.names = F)
  
  otu <- otu[, colSums(otu)>10]
  
  otu <- data.frame(otu, check.names = F)
  colors <- c("#4E79A7FF","#A0CBE8FF","#F28E2BFF","#FFBE7DFF","#59A14FFF",
              "#8CD17DFF","#B6992DFF","#F1CE63FF","#499894FF","#86BCB6FF",
              "#E15759FF","#FF9D9AFF","#79706EFF","#BAB0ACFF","#D37295FF",
              "#FABFD2FF","#B07AA1FF","#D4A6C8FF","#9D7660FF","#D7B5A6FF")
  if(display == "group") {
    dat <- profile_smp2grp(otu, group, group_colnames = group_colnames) %>% profile_transRA()
    vec <- names(sort(rowSums(dat), decreasing = T))[1:(top_n - 1)]
    dat <- rownames_to_column(dat, var = "name") %>% 
      mutate(name = ifelse(name %in% vec, name, "Other")) %>% 
      group_by(name) %>% summarise_all(sum) %>% column_to_rownames(var = "name")
    if(is.null(taxa_order)) taxa_order <- names(sort(rowSums(dat), decreasing = T))
    if(is.null(taxa_color)) taxa_color <- rep(colors, time = ceiling(nrow(dat)/20))[1:nrow(dat)]
    if(is.null(group_order)) group_order <- names(sort(unlist(dat[taxa_order[1],])))
    if(is.null(plot_title)) plot_title = "Composition analysis"
    plotdat <- rownames_to_column(dat, var = "name") %>%  gather(., key = "group", value = "val", -name) %>% 
      mutate(group = factor(group, group_order), name = factor(name, taxa_order))
    p <- ggplot(plotdat, aes(group, val, fill = name)) +
      geom_bar(stat = "identity", position = position_stack(), color = "#000000", linewidth = .2, width = .8) +
      scale_fill_manual(values = taxa_color) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "", y = "Relative Abundance (%)", title = plot_title, fill = "taxa") +
      theme_classic() +
      theme(axis.line = element_line(linewidth = .4, color = "#000000"),
            axis.ticks = element_line(linewidth = .4, color = "#000000"),
            axis.text = element_text(size = 8, color = "#000000"),
            axis.title = element_text(size = 8, color = "#000000"),
            plot.title = element_text(size = 10, color = "#000000"),
            legend.text = element_text(size = 8, color = "#000000", face = "italic"),
            legend.title = element_text(size = 10, color = "#000000"),
            panel.grid = element_blank())
    message("  ggsave(file = \"plot.compos_stacked_barplot.pdf\", width = 6, height = 4)")
  } else if(display == "sample") {
    # dat <- profile_transRA(otu)
    dat <- otu
    vec <- names(sort(rowSums(dat), decreasing = T))[1:(top_n - 1)]
    dat <- rownames_to_column(dat, var = "name") %>% 
      mutate(name = ifelse(name %in% vec, name, "Other")) %>% 
      group_by(name) %>% 
      summarise_all(sum) %>% 
      column_to_rownames(var = "name")
    if(is.null(taxa_order)) taxa_order <- names(sort(rowSums(dat), decreasing = T))
    if(is.null(taxa_color)) taxa_color <- rep(colors, time = ceiling(nrow(dat)/20))[1:nrow(dat)]
    if(is.null(sample_order)) sample_order <- names(sort(unlist(dat[taxa_order[1],])))
    if(is.null(plot_title)) plot_title = "Composition analysis"
    fill_title = paste0(plot_title, " taxa")
    plotdat <- rownames_to_column(dat, var = "name") %>% 
      gather(., key = "sample", value = "val", -name) %>% 
      mutate(sample = factor(sample, sample_order), name = factor(name, taxa_order))
    p <- ggplot(plotdat, aes(sample, val, fill = name)) +
      geom_bar(stat = "identity", position = position_stack(), color = "#000000", linewidth = .2, width = 1) +
      scale_fill_manual(values = taxa_color) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(x = "", y = "Relative Abundance (%)", title = plot_title, fill = "taxa") +
      theme_classic() +
      theme(axis.line = element_line(linewidth = .4, color = "#000000"),
            axis.ticks = element_line(linewidth = .4, color = "#000000"),
            axis.ticks.x = element_blank(),
            axis.text = element_text(size = 8, color = "#000000"),
            axis.text.x = element_text(size = 6, color = "#000000", angle = 90, hjust = 1, vjust = .5),
            axis.title = element_text(size = 8, color = "#000000"),
            plot.title = element_text(size = 10, color = "#000000"),
            legend.text = element_text(size = 8, color = "#000000", face = "italic"),
            legend.title = element_text(size = 10, color = "#000000"),
            panel.grid = element_blank())
    message("  ggsave(file = \"plot.compos_stacked_barplot.pdf\", width = 8, height = 4)")
  } else {
    stop("error parameter.")
  }
  return(p)
}



# plot_taxa_boxplot ####
# taxa的差异
# 必须对接taxa.compos.tib(display = "sample")
# group_list指定的分组的因子顺序
# group_color指定分组的填充颜色
# 调用diff.test()函数计算组间差异，可选的方法有wilcox和t检验
# 输出ggplot对象的list
plot_taxa_boxplot <- function(profile, group, group_order = NULL, group_color = NULL, group_colnames = NULL, method = "wilcox",
                              xlab = "", ylab = "Relative Abundance(%)", legend_title = "group", 
                              show_legend = T, aspect_ratio = 1){
  profile <- data.frame(profile, check.names = F)
  if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_colnames)) stop("group field (sample|group)")
  if (!is.null(group_colnames)) group <- dplyr::rename(data.frame(group, check.names = F), all_of(group_colnames))
  dat <- data.frame(t(profile), check.names = F) %>% rownames_to_column(var = "sample") %>% 
    dplyr::filter(sample %in% group$sample) %>% merge(., group, by = "sample", all.x = T)
  
  # group的order
  if(is.null(group_order)) group_order <- unique(dat$group)
  # group的配色
  if(is.null(group_color)){
    color <- c("#1f78b4","#33a02c","#e31a1c","#ff7f00","#6a3d9a","#ffff99","#b15928","#a6cee3","#b2df8a","#fb9a99","#fdbf6f","#cab2d6")
    group_color <- rep(color, times = ceiling(length(group_order)/12))[1:length(group_order)]
  }
  name <- setdiff(colnames(dat), c("sample", "group"))
  
  p_list <- list()
  for (i in name) {
    dat_i <- dat %>% select(sample, group, all_of(i)) %>% 
      rename(value = all_of(i))
      
    source("F:/code/R_func/calcu_diff.R")
    comparison <- calcu_diff(select(dat_i, -group), select(dat_i, -value), group_order = group_order, method = method) %>% 
      filter(pval < 0.05) %>% 
      pull(gp) %>% 
      strsplit(x = ., split = "_vs_")
    
    p <- ggplot(dat_i, aes(x = factor(group, group_order), y = value, color = factor(group, group_order))) + 
      geom_boxplot(fill = "transparent", outlier.size = .7, lwd = .4) +
      geom_jitter(size = .7, width = .3) +
      scale_color_manual(values = group_color) +
      labs(x = xlab, y = ylab, color = legend_title, title = i) +
      ggpubr::stat_compare_means(comparisons = comparison, method = method, size = 3,
                                 method.args = list(exact = F), label = "p.signif", 
                                 tip.length = .01, step.increase = .03, vjust = .9) +
      theme_classic() +
      theme(axis.line = element_line(linewidth = .4, color = "#000000"),
            axis.ticks = element_line(linewidth = .4, color = "#000000"),
            axis.text = element_text(size = 8, color = "#000000"),
            axis.title = element_text(size = 8, color = "#000000"),
            plot.title = element_text(size = 8, color = "#000000", hjust = .5, face = "italic"),
            legend.text = element_text(size = 8, color = "#000000"),
            legend.title = element_text(size = 8, color = "#000000"),
            panel.grid = element_blank(),
            aspect.ratio = aspect_ratio) +
      guides(color = "none")
    if (isFALSE(show_legend)) {
      p <- p + guides(color = "none")
    }
    p_list[[i]] <- p
  }
  return(p_list)
}
# lapply(p, function(x) ggsave(x, filename = paste0("family_", as.character(x$labels[4]), ".pdf"), width = 4.5, height = 4.5))
