#### info ####
# Encoding: utf-8
# Author: Jinxin Meng
# Created Data：2022-5-29
# Modified Date: 2022-9-1

library(randomForest)
library(dplyr)
library(tibble)
library(tidyr)

#### Function1 K折随机森林建模 ####
# K折随机森林建模，返回一个数据框，输出一个预测数据框
# 输入一个otu表
# group文件是样本的分组信息，分组的列名为sample和group
# k为几折检验
# seed设置随机种子, ntree设置随机森林的树数量
# make.names修改一些复杂feature的名字，主要是预防报错，二者可以降低内存。我这里就是把名字替换了，R里可以直接用make.names()函数去修改某一个表的内容
rf_Kfold <- function(otu, group, k, seed = 2023, ntree = 1000){
  start_time <- Sys.time()
  otu <- data.frame(t(data.frame(otu, check.names = F)))
  # 随机按照K折采样
  n_sample <- nrow(otu)
  size <- round(n_sample/k, 0)
  count <- seq_len(n_sample)
  sample_result <- list()
  for (i in 1:(k-1)) {
    set.seed(seed)
    sample_i <- sample(x = count, size = size, replace = F)
    sample_result[[i]] <- sample_i
    count <- setdiff(count, sample_i)
  }
  sample_result[[i+1]] <- count
  # 训练和测试
  colnames(otu) <- paste0("V", seq_len(ncol(otu)))
  otu$group <- as.factor(group$group[match(rownames(otu), group$sample)])
  pred <- rbind()
  pb <- txtProgressBar(style = 3)
  for(j in 1:length(sample_result)){ # 划分测试集和训练集
    sample_j <- sample_result[[j]]
    test <- otu[sample_j,]
    train <- otu[-sample_j,]
    set.seed(seed)
    rf_model <- randomForest(group ~ ., data = train, ntree = ntree, importance = F, proximity = T)
    pred_i <- data.frame(predict(rf_model, test, type = 'prob'))
    pred <- rbind(pred, pred_i)
    setTxtProgressBar(pb, j/length(sample_result)) # 进度计算
  }
  close(pb)
  pred <- rownames_to_column(pred, var = "sample")
  end_time <- Sys.time()
  run_time <- end_time - start_time
  cat(paste0("  Time cost: ", round(as.numeric(run_time), 4), " secs.\n"))
  return(pred)
}