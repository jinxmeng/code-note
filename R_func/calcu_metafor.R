# Jinxin Meng, 2022-9-14, 2024-5-26 -----------------------------------------------------------------

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(metafor)

# meta.metafor用于小数据的meta分析 -----------------------------------------------------------------------
# 数据转换可选用反正正弦平方根变换
# 原始数据X的平方根反正弦值作为新的分析数据。当数据偏离正态较为明显时，通过X的平方根反正弦变换可使资料接近正态分布，并达到方差齐性的要求。X值的范围是0-1。   
# 反正弦变换函数类似于logit变换或log变换，这种反向作用扩大了可变范围，同时将其向中心挤压，使极端情况更容易看到。
# otu <- asin(sqrt(otu)) # 每个样本的总丰度需为1
# https://www.programmingr.com/tutorial/arcsine-transformation/
# https://bookdown.org/robcrystalornelas/meta-analysis_of_ecological_data/#recommended-citation
# ......
# dt作为输入数据，第一列是sample，第二到n为特征列，最后两列是项目内分组和项目间分组, 如下
# sample  feat1  feat2  ..  group proj  
# S1      1.0    4.8    ..  ctr   Mengjx_2021    
# S2      5.2    3.1    ..  ctr   Zhangy_2019
# S3      30.1   16.3   ..  case  Mengjx_2021
# S4      21.4   14.4   ..  case  Zhangy_2019
# ......
# group指定group这列的列名
# group_pair指定group的分组向量
# proj指定proj这列的列名
# measure指定测量效应值选的方法
# method指定计算综合效应案例间方差的方法
# 输出为数据框
# proj  d_Mean  d_Sd  d_N c_Mean  c_Sd  c_N yi  vi  measure model method_tau2 val_tau2  I2  Q Q_pval  feature estimate  ci_lb ci_ub pval
# yi: 效应值
# vi: 案例内方差
# measure: SMD ==> Hedges'g 效应值的计算方法
# val_tau2: 案例间方差的值
# I2: 案例间差异大小占总差异的比例
# Q: # 异质性检验
# Q_pval: # 异质性检验P值 越显著异质性越大
# estimate, ci_lb, ci_ub, pval: 模型的综合评估值μ, 上下限和p值

metafor_fit <- function(dt, group = "group", group_pair = c("Disease", "Control"), proj = "proj", 
                         measure = "SMD", method = "REML") {
  # 表格处理
  dt <- dt %>% 
    rename(proj = all_of(proj), group = all_of(group))
  # feature向量
  feature <- setdiff(colnames(dt), c("sample", "group", "proj"))
  meta_outp <- rbind()
  # 循环每个feature
  for (i in feature) {
    tib <- dt %>% 
      subset(group%in%group_pair[1]) %>% 
      select(all_of(i), proj) %>% 
      rename(index = all_of(i)) %>% 
      group_by(proj) %>% 
      summarise(d_Mean = mean(index), d_Sd = sd(index), d_N = n())
    tib2 <- dt %>% 
      subset(group%in%group_pair[2]) %>% 
      select(all_of(i), proj) %>% 
      rename(index = all_of(i)) %>% 
      group_by(proj) %>% 
      summarise(c_Mean = mean(index), c_Sd = sd(index), c_N = n())
    # 合并数据  
    meta_in <- merge(tib, tib2, by= "proj")
    # Calculate effect size and variance in each project. 
    # 计算效应值和案例内方差。
    # We select the method of standardized mean difference provided by Hedges. 
    # 使用Hedges提供的SMD方法计算效应值（yi）和案例内方差（vi)。
    smd_meta <- escalc(measure = measure, data = meta_in, append = T,
                       m1i = d_Mean, m2i = c_Mean, 
                       sd1i = d_Sd, sd2i = c_Sd, 
                       n1i = d_N, n2i = c_N)
    # Calculate cumulative effect size using Random-effect model. 
    # 计算累积效应值，使用随机效应模型，随机效应模型除了随机因素引起的误差外，还考虑一些案例间的差异。
    # We calculate between-case variance using REML method (restricted maximum likelihood estimator) 
    # 使用REML方法计算案例间方差。 
    # tau^2：是案例间方差，认为不同研究之间有一些其他因素导致的差异，总异质性。
    # I^2：去判断案例间差异大小占总差异的指标之一，但是I2不可以作为选择哪种模型（固定vs.随机）的依据。
    # Qt：效应值总体的异质性，是评价效应值的差异程度，表示效应值偏离均值的程度。
    # Qt越大，则效应值越离散，暗示我们有些因素对效应值有强烈的影响，我们可以去寻找一些因素，例如年龄性别，收集数据进行下一步分析。
    # Qt的优势是可以进行显著性检验的。如果p值不显著，那么我们认为案例间的差异是随机因素造成的，这种情况下就不需要往下继续进行Meta分析。
    # estimate 累积效应值，到底是大于0还是小于0，就能知道某种处理下，feature多了还是少了，是否显著，疾病下是否对feature影响明显呢。
    non_na <- smd_meta %>% dplyr::filter(!is.na(yi)) %>% nrow()
    if(non_na != 0){
        smd_rma <- rma(yi, vi, method = method, data = smd_meta)
        # merge each data
        smd_meta <- smd_rma$data %>% 
          add_column(measure = measure, # 效应值的计算方法
                    model = "RM",  # 累积效应值计算模型
                    method_tau2 = method, # 随机效应模型估计案例内方差（Tau^2）的计算方法
                    val_tau2 = as.numeric(smd_rma$tau2), # 案例间方差的值
                    I2 = paste0(round(smd_rma$I2, digits = 2), "%"), # 案例间差异大小占总差异的比例
                    Q = smd_rma$QE, # 异质性检验
                    Q_pval = smd_rma$QEp, # 异质性检验P值 越显著异质性越大
                    feature = i,
                    estimate = as.numeric(smd_rma$beta),
                    ci_lb = smd_rma$ci.lb,
                    ci_ub = smd_rma$ci.ub,
                    pval = smd_rma$pval)
        meta_outp <- rbind(meta_outp, smd_meta)
    } else {
      message("Warning: ", i, "is not suitable.")
    }
  }
  return(meta_outp)
  cat(paste0("== estimate > 0, ==> ", group_pair[1]," ==\n== estimate < 0, ==> ",group_pair[2], " =="))
}

# metafor_fit用于小数据的meta分析 针对一种指标，需要计算一个表先 ----------------------------------------------------------

metafor_fit.1 <- function(dat, dat_colnames = NULL, measure = "SMD", method = "REML"){
  if (!all(c("name", "d_mean", "d_sd", "d_n", "c_mean", "c_sd", "c_n") %in% colnames(dat)) & is.null(dat_colnames)) {
    stop("dat field (name|d_mean|d_sd|d_n|c_mean|c_sd|c_n)
     😀 😀 😀 😀 😀 please refer to following dataset: 😀 😀 😀 😀 😀:
      ———————————————————————————————————————————————————————
      |     name     d_mean  d_sd   d_n c_mean  c_sd   c_n  |
      |     <chr>     <dbl> <dbl> <int>  <dbl> <dbl> <int>  |
      ———————————————————————————————————————————————————————
      |  1  HB-MSP    1307.  264.    16  1051.  146.     8  |
      |  2  HL-MP     1763.  186.    20  1051.  146.     8  |
      |  3  HL2-MP    1578.  140.    22  1051.  146.     8  |
      |  4  HN-YNBP   1354.  180.    19  1051.  146.     8  |
      |  5  IM-IMBP   1485.  135.    17  1051.  146.     8  |
      |  6  ....  ....  ....  ....  ....  ....  ....  ....  |
      ———————————————————————————————————————————————————————")
  }
  if (!is.null(dat_colnames)) dat <- data.frame(dat, check.names = F) %>% dplyr::rename(all_of(dat_colnames))
  
  smd_meta <- escalc(measure = measure, data = dat, append = T, 
                     m1i = d_mean, m2i = c_mean, 
                     sd1i = d_sd, sd2i = c_sd, 
                     n1i = d_n, n2i = c_n) %>% 
    na.omit %>% 
    as.data.frame
  
  if(sum(!is.na(smd_meta$yi)) != 0){
    smd_rma <- rma(yi, vi, method = method, data = smd_meta, control=list(stepadj = 0.5, maxiter = 10000))
    smd_meta <- smd_rma$data %>% 
      add_column(measure = "SMD", # 效应值的计算方法
                 model = "RM",  # 累积效应值计算模型
                 method_tau2 = "REML", # 随机效应模型估计案例内方差（Tau^2）的计算方法
                 val_tau2 = as.numeric(smd_rma$tau2), # 案例间方差的值
                 I2 = paste0(round(smd_rma$I2, digits = 2), "%"), # 案例间差异大小占总差异的比例
                 Q = smd_rma$QE, # 异质性检验
                 Q_pval = smd_rma$QEp, # 异质性检验P值 越显著异质性越大
                 estimate = as.numeric(smd_rma$beta),
                 ci_lb = smd_rma$ci.lb,
                 ci_ub = smd_rma$ci.ub,
                 pval = smd_rma$pval)
  } else {
    stop("All measure value are NA! STOP analysis further")
  }
  return(smd_meta)
}
