# Encoding: utf-8
# Author: Jinxin Meng
# Created Data：2022-5-29
# Modified Date: 2022-9-1

library(dplyr)
library(tibble)
library(tidyr)
library(randomForest)

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

#### Function2 K折rep次重复 ####
# K折随机森林建模，返回一个数据框，输出一个预测数据框
# otu输入的otu格式的表格, 但是要求行是样本，列是特征，也就是正常otu表格的转置
# group文件是样本的分组信息，分组的列名为sample和group
# k为几折检验
# rep为几次重复，最后会计算均值
# seed设置随机种子, ntree设置随机森林的树数量
# make.names修改一些复杂feature的名字，主要是预防报错，二者可以降低内存。我这里就是把名字替换了，R里可以直接用make.names()函数去修改某一个表的内容
rf_Kfold_rep <- function(otu, group, k, rep = 5, seed = 2023, ntree = 1000) {
  # rep次重复
  pred <- rbind()
  for (i in seq_len(rep)) {
    cat(paste0("  [", Sys.time(), "] Repeat no.", i, " time.\n"))
    pred_i <- rf_Kfold(otu, group, k = k, seed = seed + i -1, ntree = ntree)
    pred_i <- add_column(pred_i, rep = paste0("rep_", i))
    pred <- rbind(pred, pred_i)
  }
  pred <- pred %>% 
    select(-rep) %>% 
    group_by(sample) %>% 
    summarise_all(mean)
  return(pred)
}

#### Function3 多数据集交叉验证 ####
# 多数据集交叉验证
# 原理上，就是每一个数据集进行建模，然后用其他数据集进行验证
# 输入otu表，正常的otu表格列为样本，行为feature
# group文件是样本的分组信息，分组的列名为sample/group/dataset_name,第一列必须为"sample"，第二列必须为"group", dataset_name这列是可以指定名称的
# seed设置随机种子, ntree设置随机森林的树数量
# 指定交叉验证数据集的向量
# 返回一个交叉验证结果的数据框
rf_cross_dataset_vaildate <- function(otu, group, dataset_name = "dataset", dataset_list = NULL, seed = 2022, ntree = 1000) {
  # 提取多个数据集的list
  if(is.null(dataset_name)) {
    dataset_list <- group %>% 
      select(all_of(dataset_name)) %>%
      unlist() %>%
      as.character() %>%
      unique()
  } else {
    dataset_list = dataset_list
  }
  rf <- rbind() # VKH使用了BD的健康人样本
  for (i in dataset_list) {
    # 训练
    group_i <- group[group[[dataset_name]]%in%i,]
    sample_i <- group_i %>% 
      select(sample) %>% 
      data.frame(check.names = F) %>% 
      unlist() %>% 
      as.character()
    otu_i <- otu %>% 
      select(all_of(sample_i)) %>% 
      t() %>% 
      data.frame(check.names = F)
    colnames(otu_i) <- paste0("var_", seq(ncol(otu_i)))
    otu_i$group <- factor(group_i$group[match(rownames(otu_i), group_i$sample)])
    set.seed(seed)
    rf_model <- randomForest(group ~ ., data = otu_i, ntree = ntree, importance = F, proximity = T)
    # 预测
    for (j in setdiff(dataset_list, i)) {
      group_j <- group[group[[dataset_name]]%in%j,]
      sample_j <- group_j %>% 
        select(sample) %>% 
        data.frame(check.names = F) %>% 
        unlist() %>% 
        as.character()
      otu_j <- otu %>% 
        select(all_of(sample_j)) %>% 
        t() %>% 
        data.frame(check.names = F)
      colnames(otu_j) <- paste0("var_", seq(ncol(otu_j)))    
      otu_j$group <- factor(group_j$group[match(rownames(otu_j), group_j$sample)])
      pred <- predict(rf_model, otu_j, type = 'prob') %>% 
        data.frame() %>% 
        rownames_to_column(var = "sample") %>% 
        mutate(group = group_j$group[match(.$sample, group_j$sample)])
      roc <- roc(pred$group, pred[,2])
      rf <- rbind(rf, data.frame(train_dataset = i, test_dataset = j, seed = seed, ntree = ntree, AUC = as.numeric(str_replace(roc$auc, "Area under the curve: ", ""))))
    }
  }
  return(rf)
}

#### Function4 随机森林变量重要性排序 ####
# 随机森林计算重要变量，返回重要变量数据
# 输入的otu
# seed设置随机种子
# ntree设置随机森林树
rf_variable_rank <- function(otu, group, seed = 2023, ntree = 1000) {
  otu <- data.frame(t(data.frame(otu, check.names = F)), check.names = F)
  dat_rename <- data.frame(feature = colnames(otu), 
                           feature_rename = paste0("V", seq_len(ncol(otu))))
  colnames(otu) <- dat_rename$feature_rename
  otu$group <- factor(group$group[match(rownames(otu), group$sample)])
  # 建模
  set.seed(seed)
  rf_model <- randomForest(group ~ ., data = otu, ntree = ntree, importance = T, proximity = T)
  # 变量重要性
  importance <- importance(rf_model) %>% 
    data.frame(check.names = F) %>% 
    rownames_to_column(var = "variable") %>%
    mutate(variable = dat_rename$feature[match(variable, dat_rename$feature_rename)]) %>% 
    select(variable, MeanDecreaseAccuracy) %>%
    arrange(desc(MeanDecreaseAccuracy))
  return(importance)
}

#### Function5 随机森林留一法 ####
# leave-one-out method
# K折随机森林建模，返回一个数据框，输出一个预测数据框
# otu_t输入的otu格式的表格,t(otu)
# group文件是样本的分组信息，分组的列名为group
# k为几折检验
# seed设置随机种子
rf_loom <- function(otu, group, seed = 2022, ntree = 1000){
  set.seed(seed)
  otu <- data.frame(t(otu), check.names = F)
  colnames(otu) <- paste0("var_", seq_len(ncol(otu))) 
  otu$group <- as.factor(group$group[match(rownames(otu),group$sample)])
  # 新建进度条
  pb <- txtProgressBar(style = 3)
  star_time <- Sys.time()
  pred <- rbind()
  for(i in 1:length(rownames(otu))){
    fold_test <- otu[i,]   # 取1个样品作为测试集 
    fold_train <- otu[-i,]   # 剩下的样品作为训练集
    rf_model <- randomForest(group ~ ., data = fold_train, ntree = ntree, importance = T, proximity = TRUE)
    temp <- data.frame(predict(rf_model, fold_test, type = 'prob'))
    pred <- rbind(pred, temp) 
    setTxtProgressBar(pb, i/length(rownames(otu))) # 进度计算
  }
  pred <- pred %>% 
    rownames_to_column(var = "sample")
  end_time <- Sys.time()
  close(pb)
  run_time <- end_time - star_time
  cat(paste0("Run times: ", run_time, "\n"))
  return(pred)
}

#### Function6 otu_x建模，otu_y验证 ####
# 两个otu表均是转置后的，第一个是发现集，第二个为验证集
# group文件是样本的分组信息，分组的列至少有sample和group,第一列必须为"sample"，第二列必须为"group"
# seed设置随机种子, ntree设置随机森林的树数量
# 返回一个预测结果的数据框
rf_next_vaildate <- function(otu_x, otu_y, group, seed = 2022, ntree = 1000, label = "otu_x for modeling and otu_y for predicting"){
  # record feature information
  var <- colnames(otu_x)
  # discovery data modeling
  colnames(otu_x) <- paste0("var_", seq_len(ncol(otu_x)))
  otu_x$group <- as.factor(group$group[match(rownames(otu_x), group$sample)])
  set.seed(seed)
  rf_model <- randomForest(group ~ ., data = otu_x, ntree = ntree, importance = F, proximity = T)
  # validation data predicting
  otu_y <- otu_y[var]
  colnames(otu_y) <- paste0("var_", seq(ncol(otu_y)))    
  otu_y$group <- factor(group$group[match(rownames(otu_y), group$sample)])
  pred <- predict(rf_model, otu_y, type = 'prob') %>% 
    data.frame() %>% 
    rownames_to_column(var = "sample")
  return(pred)
}