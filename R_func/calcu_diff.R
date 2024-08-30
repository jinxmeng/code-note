# Jinxin Meng, mengjx855@163.com, 20220101, 20240714, version: 0.5 -------------

# 2022-06-01: 可选择"wilcox rank-sum","one-way anova","student's t test"三种方法做差异分析；
# 2023-01-17: diff_test_profile函数对feature进行差异分析，输入的是标准otu表和group表
# 2023-12-04: 修改diff_test_profile函数中的for循环，使用purrr::map_dfr，速度上稍微快一丢丢。
# 2023-12-04: 修改diff_test函数中的rbind()，使用tibble::add_column()，速度上稍微快一丢丢。

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

# calcu_diff -------------
# difference test；
# 用于多组数据之间两两做差异分析；
# dat输入长数据，一列sample, 一列要比较的数据；
# group输入group信息，一般是第一列是sample，第二列是group；
# group_order输入一个参与差异检验的向量, 支持多个组比较，默认时所有组；
# sample  index1
# s1      20
# s2      31 
# s3      15
# s4      12
# sample  group
# s1      g1
# s2      g1 
# s3      g2
# s4      g2
calcu_diff <- function(dat, group, dat_colnames = NULL, group_colnames = NULL,
                       group_order = NULL, center_group = NULL, padj = F,
                       method = "wilcox", paired = F, plab = F, cat = NULL) {
  if (!all(c("sample", "value") %in% colnames(dat)) & is.null(dat_colnames)) stop("group field (sample|value)")
  if (!is.null(dat_colnames)) dat <- data.frame(dat, check.names = F) %>% dplyr::rename(all_of(dat_colnames))
  if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_colnames)) stop("group field (sample|group)")
  if (!is.null(group_colnames)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_colnames))
  if (is.null(group_order)) group_order <- unique(group$group)
  if (!method %in% c("wilcox", "anova", "t")) stop("method option: wilcox, anova, t")
  if (!is.null(center_group)) {
    if (!(center_group %in% group_order)) {
      stop("center_group not existing.")
    }
  } 
  
  dat <- mutate(dat, group = group$group[match(sample, group$sample)])
  
  if (!is.null(cat)) cat(paste0("    Processing: ", cat, "\n"))
  pb <- txtProgressBar(style = 3, width = 50, char = "#")
  diff <- data.frame(gp = character(), pval = numeric(), method = character())
  if (!is.null(center_group)) {
    group_order <- setdiff(group_order, center_group)
    if (method == "wilcox") {
      dat_center <- subset(dat, group %in% center_group) %>% dplyr::select(value) %>% unlist() %>% as.numeric()
      for (i in 1:length(group_order)) {
        dat_i <- subset(dat, group %in% group_order[i]) %>% dplyr::select(value) %>% unlist() %>% as.numeric()
        wilcox <- wilcox.test(dat_i, dat_center, paired = paired, exact = F)
        gp <- paste0(group_order[i], "_vs_", center_group)
        diff <- tibble::add_row(diff, gp = gp, pval = as.numeric(wilcox$p.value), method = "wilcoxon-rank sum")
        setTxtProgressBar(pb, i/length(group_order))
      }
    } else if (method == "anova") {
      for (i in 1:length(group_order)) {
        group_i <- subset(group, group %in% c(group_order[i], center_group))
        dat_i <- subset(dat, sample %in% group_i$sample)
        anova <- oneway.test(value ~ group, dat_i)
        gp <- paste0(group_order[i], "_vs_", center_group)
        diff <- tibble::add_row(diff, gp = gp, pval = as.numeric(anova$p.value), method = "oneway anova")
        setTxtProgressBar(pb, i/length(group_order))
      }
    } else if (method == "t") {
      for (i in 1:length(group_order)) {
        group_i <- subset(group, group %in% c(group_order[i], center_group))
        dat_i <- subset(dat, sample %in% group_i$sample)
        t <- stats::t.test(value ~ group, dat_i, paired = paired)
        gp <- paste0(group_order[i], "_vs_", center_group)
        diff <- tibble::add_row(diff, gp = gp, pval = as.numeric(t$p.value), method = "student's t")
        setTxtProgressBar(pb, i/length(group_order))
      }
    }
  } else {
    x = 0
    if (method == "wilcox") {
      for (i in 1:(length(group_order) - 1)) {
        for (j in (i + 1):length(group_order)) {
          x = x + 1
          dat_i <- subset(dat, group %in% group_order[i]) %>% dplyr::select(value) %>% unlist() %>% as.numeric()
          dat_j <- dat %>% subset(group %in% group_order[j]) %>% dplyr::select(value) %>% unlist() %>% as.numeric()
          wilcox <- wilcox.test(dat_i, dat_j, paired = paired, exact = F)
          gp <- paste0(group_order[i], "_vs_", group_order[j])
          diff <- tibble::add_row(diff, gp = gp, pval = as.numeric(wilcox$p.value), method = "wilcoxon-rank sum")
          setTxtProgressBar(pb, x/choose(length(group_order), 2))
        }
      }
    } else if (method == "anova") {
      for (i in 1:(length(group_order) - 1)) {
        for (j in (i + 1):length(group_order)) {
          x = x + 1
          group_ij <- subset(group, group %in% c(group_order[i], group_order[j]))
          dat_ij <- subset(dat, sample %in% group_ij$sample)
          anova <- oneway.test(value ~ group, dat_ij)
          gp <- paste0(group_order[i], "_vs_", group_order[j])
          diff <- tibble::add_row(diff, gp = gp, pval = as.numeric(anova$p.value), method = "oneway anova")
          setTxtProgressBar(pb, x/choose(length(group_order), 2))
        }        
      } 
    } else if (method == "t") {
      for (i in 1:(length(group_order) - 1)) {
        for (j in (i + 1):length(group_order)) {
          x = x + 1
          group_ij <- subset(group, group %in% c(group_order[i], group_order[j]))
          dat_ij <- subset(dat, sample %in% group_ij$sample)
          t <- stats::t.test(value ~ group, dat_ij, paired = paired)
          gp <- paste0(group_order[i], "_vs_", group_order[j])
          diff <- tibble::add_row(diff, gp = gp, pval = as.numeric(t$p.value), method = "student's t")
          setTxtProgressBar(pb, x/choose(length(group_order), 2))
        }
      }
    }
  }
  if (isTRUE(plab)) diff <- add_column(diff, plab = ifelse(diff$pval < 0.001, "***", ifelse(diff$pval < 0.01, "**", ifelse(diff$pval < 0.05, "*", ""))), .after = "pval")
  if (isTRUE(padj)) diff <- add_column(diff, padj = p.adjust(diff$pval, method = "BH"), .after = "pval")
  close(pb)
  return(diff)
}

# calcu_diff_profile ---------------
# 用于对OTU表中所有的feature进行差异分析；
# profile输入OTU表， 行为feature，列为sample；
# group输入group表，一般是第一列是sample，第二列是group；如果不是这种情况，请用group_colnames指定一个改名字的向量
# group_order输入一个参与差异检验的向量, 支持多个组比较，默认时所有组
#           s1  s2  s3  s4
# feature1  20  31  15  12
# feature2  21  32  16  13
# feature3  22  33  17  14
# feature4  23  34  18  15
# sample  group
# s1      g1
# s2      g1 
# s3      g2
# s4      g2
calcu_diff_profile <- function(profile, group, group_order = NULL, group_colnames = NULL, 
                               center_group = NULL, method = "wilcox", paired = F, 
                               plab = F, padj = F) {
  if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_colnames)) stop("group field (sample|group)")
  if (!is.null(group_colnames)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_colnames))
  if (is.null(group_order)) group_order <- unique(group$group)
  
  dat <- data.frame(t(profile), check.names = F) %>% 
    rownames_to_column(var = "sample")
  tmp_vec <- setdiff(colnames(dat), "sample")
  diff <- purrr::map_dfr(tmp_vec, \(x) 
                         calcu_diff(dat = dplyr::select(dat, sample, value = all_of(x)), 
                                    group = group, group_order = group_order, method = method, 
                                    paired = paired, center_group = center_group, cat = x) %>% 
                           add_column(name = x, .before = "gp")
                         )
  if (isTRUE(padj)) diff <- add_column(diff, padj = p.adjust(diff$pval, method = "BH"), .after = "pval")
  if (isFALSE(padj) & isTRUE(plab)) diff <- add_column(diff, plab = ifelse(diff$pval < 0.001, "***", ifelse(diff$pval < 0.01, "**", ifelse(diff$pval < 0.05, "*", ""))), .after = "pval")
  if (isTRUE(padj) & isTRUE(plab)) diff <- add_column(diff, plab = ifelse(diff$padj < 0.001, "***", ifelse(diff$padj < 0.01, "**", ifelse(diff$padj < 0.05, "*", ""))), .after = "padj")
  return(diff)
  }
