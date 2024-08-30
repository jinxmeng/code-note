#### info ####
# encoding: utf-8
# author: Jinxin Meng
# e-mail: mengjx855@163.com
# created data：2024-03-07
# modified data: 2024-03-07
# version: 0.1

library(dplyr)
library(Hmisc) %>% suppressMessages()
library(minpack.lm) %>% suppressMessages()
library(stats4) %>% suppressMessages()

# calcu_NCM ----
calcu_NCM <- function(profile, name = NULL) {
  # 将 Sloan 等（2006）的中性模型拟合到一个物种或分类群的丰度表，并返回几个拟合统计数据。
  # 或者，将根据它们在元群落中的丰度返回每个分类群的预测出现频率
  # 用非线性最小二乘法（Non-linear least squares，NLS）拟合模型参数
  spp <- t(profile)
  N <- mean(apply(spp, 1, sum)) # 计算总相对丰度的平均值
  p.m <- apply(spp, 2, mean) # # 计算每个物种的平均相对丰度
  p.m <- p.m[p.m != 0] # 去除平均值为0的物种
  p <- p.m / N  # 计算每个物种的相对丰度
  spp.bi <- 1 * (spp > 0) # 将原始数据二值化，表示物种的存在与否
  freq <- apply(spp.bi, 2, mean) # # 计算每个物种的出现频率
  freq <- freq[freq != 0] # 去除频率为0的物种
  C <- merge(p, freq, by=0) # 合并相对丰度和频率数据
  C <- C[order(C[,2]),] # 按照频率排序
  C <- as.data.frame(C) # 将结果转为数据框
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] # 去除包含0的行
  p <- C.0[,2] # 提取相对丰度和频率的数据
  freq <- C.0[,3] # 提取相对丰度和频率的数据
  names(p) <- C.0[,1] # 为数据命名
  names(freq) <- C.0[,1] # 为数据命名
  d = 1 / N # 计算d的值
  m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = F),start = list(m = 0.1)) # 使用非线性最小二乘法（NLS）拟合模型参数 m（或 Nm）
  m.fit  #获取 m 值
  m.ci <- confint(m.fit, 'm', level = 0.95) # 计算 m 的置信区间
  freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = F) # 计算预测的频率值
  pred.ci <- binconf(freq.pred * nrow(spp), nrow(spp), alpha = 0.05, method = "wilson", return.df = T) # 使用 Wilson 方法计算置信区间
  Rsqr <- 1 - (sum((freq - freq.pred) ^ 2)) / (sum((freq - mean(freq)) ^ 2)) # 计算 R2 值
  res <- data.frame(Rsqr = Rsqr, m = as.numeric(coef(m.fit)), mN = as.numeric(coef(m.fit) * N), model = "NCM")
  if (!is.null(name)) res$name = name
  return(res)
}

# plot_NCM ----
# ref: 棱镜TuT & 小白鱼的生统笔记
# input a profile
plot_NCM <- function(profile, title = NULL) {
  # 将 Sloan 等（2006）的中性模型拟合到一个物种或分类群的丰度表，并返回几个拟合统计数据。
  # 或者，将根据它们在元群落中的丰度返回每个分类群的预测出现频率
  # 用非线性最小二乘法（Non-linear least squares，NLS）拟合模型参数
  spp <- t(profile)
  N <- mean(apply(spp, 1, sum)) # 计算总相对丰度的平均值
  p.m <- apply(spp, 2, mean) # # 计算每个物种的平均相对丰度
  p.m <- p.m[p.m != 0] # 去除平均值为0的物种
  p <- p.m / N  # 计算每个物种的相对丰度
  spp.bi <- 1 * (spp > 0) # 将原始数据二值化，表示物种的存在与否
  freq <- apply(spp.bi, 2, mean) # # 计算每个物种的出现频率
  freq <- freq[freq != 0] # 去除频率为0的物种
  C <- merge(p, freq, by=0) # 合并相对丰度和频率数据
  C <- C[order(C[,2]),] # 按照频率排序
  C <- as.data.frame(C) # 将结果转为数据框
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] # 去除包含0的行
  p <- C.0[,2] # 提取相对丰度和频率的数据
  freq <- C.0[,3] # 提取相对丰度和频率的数据
  names(p) <- C.0[,1] # 为数据命名
  names(freq) <- C.0[,1] # 为数据命名
  d = 1 / N # 计算d的值
  m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = F),start = list(m = 0.1)) # 使用非线性最小二乘法（NLS）拟合模型参数 m（或 Nm）
  m.fit  #获取 m 值
  m.ci <- confint(m.fit, 'm', level = 0.95) # 计算 m 的置信区间
  freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = F) # 计算预测的频率值
  pred.ci <- binconf(freq.pred * nrow(spp), nrow(spp), alpha = 0.05, method = "wilson", return.df = T) # 使用 Wilson 方法计算置信区间
  Rsqr <- 1 - (sum((freq - freq.pred) ^ 2)) / (sum((freq - mean(freq)) ^ 2)) # 计算 R2 值
  Rsqr  #获取模型的 R2
  
  # 输出 3 个统计结果数据表，包括各物种或分类群的平均相对丰度（p.csv）、出现频率（freq.csv）和预测的出现频率（freq.pred.csv）
  # write.csv(p, file = "p.csv") # p 是平均相对丰度（mean relative abundance）
  # write.csv(freq, file = "freq.csv") # freq 是出现频率（occurrence frequency）的观测值
  # write.csv(freq.pred, file = "freq.pred.csv") # freq.pred 是出现频率（occurrence frequency）的预测值，即中性模型的拟合值
  
  plotdat <- data.frame(p, freq, freq.pred, pred.ci[,2:3]) %>% 
    dplyr::mutate(class = ifelse(freq <= Lower, "lower", ifelse(freq >= Upper, "upper", "inter")))
  colors <- structure(c("#000000", "#a52a2a", "#29a6a6"), names = c("inter", "lower", "upper"))
  subtitle <- paste0("Coefficient: Rsqr = ", round(Rsqr, 4), ", m = ", round(coef(m.fit), 4), ", mN = ", round(coef(m.fit) * N, 4))
  if (is.null(title)) title <- "Neutral community model analysis"

  p <- ggplot(plotdat) +
    geom_point(aes(x = log10(p), y = freq, color = class), size = 1.2, show.legend = F) +
    geom_line(aes(x = log10(p), y = freq.pred, group = 1), linewidth = .6, color = "blue") +
    geom_line(aes(x = log10(p), y = Upper, group = 1), linewidth = .6, color = "grey60", lty = "dashed") +
    geom_line(aes(x = log10(p), y = Lower, group = 1), linewidth = .6, color = "grey60", lty = "dashed") +
    scale_color_manual(values = colors) +
    labs(x = "Mean Relative Abundance (log10)", y = "Frequency of Occurance",
         title = title, subtitle = subtitle) +
    theme_bw() +
    theme(axis.line = element_blank(),
          axis.ticks = element_line(linewidth = .4, color = "#000000"),
          axis.text = element_text(size = 8, color = "#000000"),
          axis.title = element_text(size = 8, color = "#000000"),
          plot.title = element_text(size = 10, color = "#000000"),
          plot.subtitle = element_text(size = 8, color = "#000000"),
          legend.text = element_text(size = 8, color = "#000000"),
          legend.title = element_text(size = 8, color = "#000000"),
          panel.grid = element_blank(),
          aspect.ratio = 3/4)
  message("  ggsave(file = \"plot.NCM.pdf\", width = 4.5, height = 4)")
  return(p)
}