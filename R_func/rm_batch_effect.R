# Encoding: utf-8
# Author: Jinxin Meng
# Email: mengjx855@163.com
# Created Data：2022-11-2
# Modified Data: 2022-11-2
# Version: 0.1

library(dplyr)
library(tibble)
library(tidyr)
library(sva)
library(limma)

#### Function1 ####
# remove batch effect
# ComBat() from sva package
# 目前只支持一个实验分组的去批次效应
# profile输入一个expr.profile或OTU.profile
# design输入一个分组表要带有一个批次列
# sample  group batch
# s1    case    bat_1
# s2    ctr     bat_1
# s3    case    bat_2
# s4    ctr     bat_2
# BiocManager::install("sva", force = T)
rm_batch_combat <- function(profile, design, design_n = "group", batch_n = "batch", par.prior = T) {
  # 数据整理
  design <- design[match(colnames(profile), design$sample),] %>%
    data.frame(row.names = NULL) %>% 
    rename(group = all_of(design_n),
           batch = all_of(batch_n)) %>% 
    mutate(group = factor(group),
           batch = factor(batch)) %>%
    column_to_rownames(., var = "sample")
  # Sva包消除批次效应
  # 建立批次效应的模型，design表示的是数据中除了有不同的批次，还有生物学上的差异。
  # 在校正的时候要保留生物学上的差异，不能矫正过枉。
  model <- model.matrix(~ as.factor(group), data = design)
  # 消除批次效应
  # combat_Expr就是校正后的数据
  # par.prior T - 参数检验 F - 非参数检验，初步理解为根据正态分布选择
  profile2 <- ComBat(dat = profile, batch = design$batch, mod = model, par.prior = par.prior)
  return(profile2)
}

#### Function2 ####
# remove batch effect
# removeBatchEffect() from limma package
# 目前只支持一个实验分组的去批次效应
# profile输入一个expr.profile或OTU.profile
# design输入一个分组表要带有一个批次列，多了几列没关系，设置好batch_n就行
# sample  batch
# s1  bat_1
# s2  bat_1
# s3  bat_2
# s4  bat_2
rm_batch_effect <- function(profile, design, batch_n = "batch") {
  # 数据整理
  batch <- design[match(colnames(profile), design$sample),] %>% 
    rename(batch = all_of(batch_n)) %>% 
    select(batch) %>% 
    unlist() %>% 
    as.character() 
  # 去除批次效应
  profile2 <- removeBatchEffect(profile, batch)
  return(profile2)
}