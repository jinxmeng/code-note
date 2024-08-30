# Encoding: utf-8
# Author: Jinxin Meng
# Email: mengjx855@163.com
# Created Data：2023-11-10
# Modified Data: 2023-12-16
# Version: 0.1

pacman::p_load(dplyr, tibble, tidyr)

#### Function 1 ####
# 处理psych corr.test 分析结果
# 保留|cor|≥rho且pval<pval的值
# 若一行/列中没有符合条件的相关性值，删除这行/列.
# 可以输出data.frame格式，会进一步过滤结果
# 可以全部输出
corr_process <- function(psych_corr_object, rho = 0.5, pval = NULL, padj = NULL, out_df = F, out_all = F) {
  if (isTRUE(out_all)) {
    pval = 1
    padj = 1
    rho = 0
  }
  corr <- psych_corr_object
  corr_r <- as.matrix(corr$r)
  if (is.null(pval) & is.null(padj)) {
    p_cutoff <- 0.05
    p_title <- "padj"
    corr_p <- as.matrix(corr$p.adj)
  } else if (is.null(pval) & !is.null(padj)) {
    p_cutoff <- padj
    p_title <- "padj"
    corr_p <- as.matrix(corr$p.adj)
  } else if (!is.null(pval) & is.null(padj)) {
    p_cutoff <- pval
    p_title <- "pval"
    corr_p <- as.matrix(corr$p)
  } else if (!is.null(pval) & !is.null(padj)) {
    p_cutoff <- padj
    p_title <- "padj"
    corr_p <- as.matrix(corr$p.adj)
  }
  tmp_r <- corr_r
  tmp_p <- corr_p
  tmp_r[abs(tmp_r) < rho] <- 0
  tmp_p[tmp_p >= p_cutoff] <- -1
  tmp_p[tmp_p < p_cutoff & tmp_p >= 0] <- 1
  tmp_p[tmp_p == -1] <- 0
  tmp_corr <- tmp_r*tmp_p
  tmp_corr <- tmp_corr[apply(tmp_corr, 1, \(x) sum(x!=0)>0), apply(tmp_corr, 2, \(x) sum(x!=0)>0)] %>% 
    data.frame(check.names = F)
  
  if (isTRUE(out_df)) {
    tmp_r <- data.frame(corr_r[rownames(tmp_corr), colnames(tmp_corr)], check.names = F) %>% 
      rownames_to_column(var = "name_x") %>% 
      gather(key = "name_y", value = "r", -name_x)
    tmp_p <- data.frame(corr_p[rownames(tmp_corr), colnames(tmp_corr)], check.names = F) %>% 
      rownames_to_column(var = "name_x") %>% 
      gather(key = "name_y", value = "pval", -name_x)
    tmp_corr <- merge(tmp_r, tmp_p, by = c("name_x", "name_y")) %>% 
      dplyr::filter(abs(r)>rho & pval < p_cutoff)
    if (p_title == "padj") {
      colnames(tmp_corr)[4] <- "padj"
    }
  }
  return(tmp_corr)
}

# corr_merge -----------------------
# merge cor pval and rho matrix.
# 保留|cor|≥rho且pval<pval的值
# 若一行/列中没有符合条件的相关性值，删除这行/列.
corr_merge <- function(r_mat, p_mat, rho = 0.5, pval = 0.05) { 
  corr_r <- as.matrix(r_mat)
  corr_p <- as.matrix(p_mat)
  corr_r[abs(corr_r) < rho] <- 0
  corr_p[corr_p >= pval] <- -1
  corr_p[corr_p < pval & corr_p >= 0] <- 1
  corr_p[corr_p == -1] <- 0
  corr <- corr_r*corr_p
  corr <- corr[apply(corr, 1, \(x) sum(x!=0)>=1), apply(corr, 2, \(x) sum(x!=0)>=1)] %>% 
    data.frame(check.names = F)
  return(corr)
}
