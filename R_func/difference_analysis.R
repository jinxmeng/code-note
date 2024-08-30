# Jinxin Meng, 20231101, 20240630 -----------

# Version: 1.3
# 2023-10-25: add xlim parameter in plot_roc.
# 2023-10-25: add trans parameter in get_feature_diff.
# 2023-10-25: add map_name parameter in the get_feature_diff.
# 2023-10-25: add names in the plot_volcano.
# 2023-11-01: add prevalence statistics in output of get_feature_diff, min_abundance as a threshold of presence.
# 2023-11-02: modify group_color in plot_volcano using function structure

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
source("F:/code/R_func/utilities.R")
source("F:/code/R_func/profile_process.R")


# difference_analysis -------------------------------------
# 对feature进行组间富集分析，绘制火山图并保存数据
# profile输入正常的profile表
# group表中分组列必须为group
# gp传入一个向量；例如c("A","B")，之后的分析即为AvsB，上调在A中富集，下调在B中富集，初步版本只能是"Disease", "Control"
# log2fc为认定富集的差异倍数阈值
# padj为认定富集的显著性阈值;
# group[data.frame]: metadata contain sample and group column, also specify by map_names parameter.
# trans[log]: used in difference analysis to transfrom profile, only including LOG in current version. [NULL, RA, LOG]
# min_abundace[num]: the threshold of a feature presencing when evaluate prevalence rate (%).

difference_analysis <- function(profile, group, gp = NULL, trans = NULL, min_abundance = 0,
                                group_colnames = NULL, diff_method = "wilcox", add_diff_lab = F, 
                                log2fc = NULL, pvalue = NULL, qvalue = NULL) {
  
  if (!all(colnames(group) %in% c("sample", "group")) & is.null(group_colnames)) stop("group field (sample|group)")
  if (!is.null(group_colnames)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_colnames))

  # difference analysis order is define.
  if (is.null(gp)) { gp <- unique(group$group)
  } else { gp <- gp }
  
  if (!all(group$group %in% gp)) {
    group <- dplyr::filter(group, group %in% gp)
    profile <- dplyr::select(profile, all_of(group$sample))
  }
    
  profile <- data.frame(profile, check.names = F)
    
  if (is.null(trans)) { 
    data <- data.frame(t(profile), check.names = F) 
  } else if (trans == "RA") { 
    data <- data.frame(t(profile_transRA(profile)), check.names = F) 
  } else if (trans == "LOG") { 
    data <- data.frame(t(profile_transLOG(profile)), check.names = F) 
  } else { 
    stop("in profile transformation, parameter only including RA, LOG in current version. 😀 ")
  }
  data <- dplyr::mutate(data, group = group$group[match(rownames(data), group$sample)])
  
  # difference analysis using diff_method.
  cat(paste0("[", Sys.time(), "] Hypothesis testing using ", diff_method, " method.\n"))
  
  name <- setdiff(colnames(data), "group")
  pb <- txtProgressBar(style = 3, width = 50, char = "#")
  pval_vec <- c()
  
  for (i in 1:length(name)) {
    data_i <- dplyr::select(data, all_of(name[i]), group)
    if (diff_method == "wilcox") {
      test <- wilcox.test(unlist(subset(data_i, group%in%gp[1])[1]),
                          unlist(subset(data_i, group%in%gp[2])[1]),
                          paired = F, exact = F)
    } else if (diff_method == "t") { # t
      test <- stats::t.test(unlist(subset(data_i, group%in%gp[1])[1]),
                            unlist(subset(data_i, group%in%gp[2])[1]), 
                            paired = F)
    }
    
    pval_vec[i] <- test$p.value
    setTxtProgressBar(pb, i/length(name)) 
  }
  close(pb)
  
  diff <- data.frame(name = name, pval = pval_vec) %>% 
    mutate(padj = p.adjust(pval, method = "BH"))
  
  # calu feature abundance and foldchange
  cat(paste0("[", Sys.time(), "] Abundance analysis.\n"))
  
  abund <- data %>%
    group_by(group) %>%
    summarise_all(mean) %>%
    column_to_rownames(var = "group") %>% 
    t(.) %>% 
    data.frame(check.names = F) %>% 
    rownames_to_column(var = "name") %>% 
    dplyr::select(name, index1_ab = all_of(gp[1]), index2_ab = all_of(gp[2])) %>% 
    mutate(comparsion = paste(gp, collapse = "_vs_"),
           FC = index1_ab / index2_ab,
           log2FC = log2(FC))
  
  # calu feature prevalence
  cat(paste0("[", Sys.time(), "] Prevalence statistic.\n"))
  
  prevalence <- data %>% 
    group_by(group) %>% 
    group_modify(~ purrr::map_df(.x,\(x) sum(x > min_abundance)/length(x) *100)) %>% 
    column_to_rownames(var = "group") %>% 
    t() %>% 
    data.frame(check.names = F) %>% 
    rownames_to_column(var = "name") %>% 
    dplyr::select(name, index1_pvl = all_of(gp[1]), index2_pvl = all_of(gp[2]))
  
  # reshape the output
  cat(paste0("[", Sys.time(), "] Reshape result.\n"))
  
  out <- merge(abund, prevalence, by = "name", all = T) %>% 
    merge(., diff, by = "name", all = T) %>% 
    relocate(index1_pvl, .after = index1_ab) %>% 
    relocate(index2_pvl, .after = index2_ab)
  
  out <- out %>% 
    dplyr::filter(!is.nan(FC)) %>% 
    within(., {
           FC[FC == Inf] <- max(FC[!is.infinite(FC)])
           log2FC[log2FC == Inf] <- max(log2FC[!is.infinite(log2FC)])
           log2FC[log2FC == -Inf] <- min(log2FC[!is.infinite(log2FC)])
           })
  
  if(isTRUE(add_diff_lab)) {
    if(is.null(log2fc)) log2fc = 1; cat(paste0("[", Sys.time(), "] Default setting: log2fc = 1\n"))
    if(is.null(pvalue) & is.null(qvalue)) pvalue = 0.05; cat(paste0("[", Sys.time(), "] Default setting: pval = 0.05\n"))
    if(!is.null(pvalue)) out <- mutate(out, enriched = ifelse(log2FC > log2fc & pval < pvalue, gp[1], ifelse(log2FC < -log2fc & pval < pvalue, gp[2], "none")))
    if(!is.null(qvalue)) out <- mutate(out, enriched = ifelse(log2FC > log2fc & padj < qvalue, gp[1], ifelse(log2FC < -log2fc & padj < qvalue, gp[2], "none")))
  }
  
  # rename the output
  names_vec <- c("index1_ab", "index2_ab", "index1_pvl", "index2_pvl")
  names(names_vec) <- c(paste0(gp[1], "_ab"), paste0(gp[2], "_ab"), 
                        paste0(gp[1], "_pvl"), paste0(gp[2], "_pvl"))
  out <- dplyr::rename(out, all_of(names_vec))
  
  cat(paste0("[", Sys.time(), "] Program end.\n"))
  return(out)
}

# add_plab -----------------------------------------
add_plab <- function(data, by = "pval", format = 1){
  if (sum(colnames(data) %in% by) == 0) stop("data must have columns (pval|padj)")
  if (format == 1) {
    data$plab <- cut(dplyr::pull(data, all_of(by)),
                     breaks = c(0, 0.001, 0.01, 0.05, 1),
                     labels = c("***", "**", "*", "")) %>% 
      as.character()
  } else if (format == 2) {
    data$plab <- cut(dplyr::pull(data, all_of(by)),
                     breaks = c(0, 0.001, 0.01, 0.05, 1),
                     labels = c("#", "+", "*", "")) %>% 
      as.character()
  } else if (format == 3) {
    data$plab <- cut(dplyr::pull(data, all_of(by)),
                     breaks = c(0, 0.001, 0.01, 0.05, 1),
                     labels = c("+", "**", "*", "")) %>% 
      as.character()
  }  
  data <- relocate(data, .after = all_of(by))
  return(data)
}

