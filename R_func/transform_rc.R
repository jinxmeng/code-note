# Encoding: utf-8
# Creator: Jinxin Meng
# Created Data：2022-05-29
# Modified Data: 2023-09-20
# Version: 1.1
# 2023-09-20: 修改了一些变量名称，增加count2fpkm的函数

library(dplyr)
library(tibble)
library(tidyr)

# file format
# profile|profile
#         s1  s2 ...
# gene1   2   3  
# gene2   12  2  
# ...
# len|length
# gene1   124
# gene2   398
# ...
# TPM的计算方法：把比对到的某个基因的fragment数目，乘以1e6，除以基因的长度，再求相对丰度；
# FPKM的计算方法：把比对到的某个基因的fragment数目，乘以1e9，除以基因的长度，其比值再除以总map的reads数；
# RPM的计算方法：把比对到的某个基因的fragment数目，乘以1e6，除以此样本fragment总数；
rc2tpm <- function(profile, len, len_colnames = NULL) {
  if (!all(colnames(len) %in% c("name", "len")) & is.null(len_colnames)) stop("len field (name|len)")
  if (!is.null(len_colnames)) dat <- data.frame(len, check.names = F) %>% dplyr::rename(all_of(len_colnames))
  len <- tibble::column_to_rownames(len, var = "name")
  len <- len[rownames(profile),]
  tpm <- apply(profile, 2, function(x) x/len) %>% 
    apply(., 2, function(x) x/sum(x)*1e6) %>% 
    as.data.frame(., check.names = F)
  return(tpm)
}


rc2fpkm <- function(profile, len, len_colnames = NULL) {
  if (!all(colnames(len) %in% c("name", "len")) & is.null(len_colnames)) stop("len field (name|len)")
  if (!is.null(len_colnames)) len <- data.frame(len, check.names = F) %>% dplyr::rename(all_of(len_colnames))
  len <- tibble::column_to_rownames(len, var = "name")
  len <- len[rownames(profile),]
  reads <- colSums(profile)
  fpkm <- apply(profile, 2, function(x) x*1e12/len) %>% 
    apply(., 1, function(x) x/reads) %>% 
    t() %>% 
    data.frame(check.names = F)
  return(fpkm)
} 