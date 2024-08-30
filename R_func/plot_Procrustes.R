#### Info ####
# Jinxin Meng, 20231215, 20240126
# Version: 1.0
# procrustes分析|普鲁克分析|两个矩阵的相关性分析
# Procrustes分析与Mantel Test的比较
# 提到群落分析，Procrustes分析与Mantel test都是用于分析物种组成和环境属性关系的常见方法，当然二者的具体关注点还是有区别的，方法各有自身的优点。
# 但如果只聚焦在评估两数据集一致性上，似乎Procrustes分析更直观一些。
# M2统计量及其显著性检验p值提供了两个数据集之间一致性的总体度量，同时数据集的图形匹配和相关残差提供了比Mantel test更丰富的信息源。
# 在对应点的坐标匹配度较好时，两个数据集表现出良好的一致性。坐标匹配度越差表明这些点与整体趋势不匹配，
# 这类似于回归分析中残差较大的点，这些点不符合样本的总体趋势。
# 此外，PROTEST的统计功效也被证明优于Mantel test的统计功效（Peres-Neto and Jackson, 2001）。
# 因此，如果两组数据之间存在潜在关系，则Procrustes分析更有能力检测到，并且鉴于结果的图形性质，它还提供了出色的解释性准则。

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(vegan)

#### plot_Procrustes ####
# Procrustes analysis 
# input: profile[feature, sample] | dist 
# 输入两个profile表，或者是基于两个profile表分析的聚类矩阵。
# 默认是y旋转变换到x
plot_Procrustes <- function(profile_x, profile_y, dist_x, dist_y, dis_method = "bray") {
  if (!missing(profile_x) & !missing(profile_y)) {
    dist_x <- vegdist(t(profile_x), method = dis_method)
    dist_y <- vegdist(t(profile_y), method = dis_method)
  } else if (missing(dist_x) & missing(dist_y)) {
    stop("distance matrix missing!")
  }
  # 降维
  PCoA_x <- cmdscale(dist_x)
  PCoA_y <- cmdscale(dist_y)
  
  # 选择模式
  # 通过对该函数中各参数的了解，可知X为目标矩阵也就是降维后的环境（功能基因等）坐标，Y为降维后的物种数据的坐标，
  # 因为后续普氏分析中旋转和缩放操作是针对Y，将Y匹配给X。
  # 另外，当symmetric=FALSE时，处于"非对称"模式，X和Y的分配值调换后，普氏分析的偏差平方和（M^2^）也会随之改变。
  # 而当symmetric=TRUE时，从而给出更合适比例的对称统计。
  # “对称”模式下，X和Y的分配值调换后，普氏分析的偏差平方和（M2）不会发生改变，但注意旋转仍将是非对称的。
  # 以对称模式进行普氏分析（symmetric = TRUE）
  # 作者回应 https://stats.stackexchange.com/questions/563911/what-is-the-difference-between-symmetric-and-non-symmetric-in-procrustes-protest
  # 将矩阵B拟合到目标矩阵A：我想获得与A类似的旋转/缩放/平移矩阵B。如果这是您的目标，您应该使用非对称(B到A)旋转。 
  # 如果您不想获得结果（旋转），但您只是对两个矩阵之间的相似性以及该相似性的统计量感兴趣，则应该使用对称旋转。 
  proc <- procrustes(PCoA_x, PCoA_y, symmetric = T)
  summary(proc)
  # # 评价一致性
  # 如果样本中物种与环境一致性（相似性）越近，则对应的残差越小，
  # 反之物种与环境的相似性越远，则残差越大（三条辅助线对应的位置分别为残差25%、50%和75%）
  # plot(proc, kind = 2)
  # residuals(proc)
  
  # 置换检验999次
  # 普氏分析中M2统计量的显著性检验
  # 在 Procrustes 分析中，M2 常常指代的是 Procrustes 统计量，它是度量两个形状间差异程度的一个指标。
  # 更具体地说，M2 是源数据集的点经过旋转、缩放和/或平移后与目标数据集中对应点的平方距离和。
  # 它代表了变换后源数据集的点与目标数据集点之间的不匹配程度。M2 的数值越小，表明两组数据集的形状越相似；反之，则表明它们之间的差异越大。
  set.seed(2024)
  proc_test <- protest(PCoA_x, PCoA_y, permutations = 999)
  # Procrustes Sum of Squares (m12 squared): M2统计值
  # Significance:  0.001  # p value
  
  # 提取结果
  # 偏差平方和（M2统计量）
  m2 = round(proc_test$ss, 4)
  # 对应p值结果
  pval = proc_test$signif
  # label
  label <- paste0("'coefficients:'~'M'^2~'='~'", m2, "'~~italic('p')~'<'~'", pval, "'") %>% as.formula() %>% eval()
  
  # 首先提取降维后的数据轴1和2的坐标，并且提取转换的坐标；然后进行绘制。
  # 获得x和y轴的坐标及旋转过的坐标
  proc_point <- cbind(
    data.frame(proc$Yrot) %>% rename(X1_rotated = X1, X2_rotated = X2),# Y-矩阵旋转后的坐标，就是物种矩阵为了逼近X做了调整，调整后的坐标。
    data.frame(proc$X) %>% rename(X1_target = Dim1, X2_target = Dim2) # X-目标矩阵，就是proc分析输入的第一个矩阵，作为目标矩阵没变化。
  )
  proc_coord <- data.frame(proc$rotation) # Y旋转后的坐标轴
  
  # 绘图
  p <- ggplot(proc_point) + # 旋转坐标到目的坐标的一半上一个颜色，另一半上不同的颜色
    geom_segment(aes(x = X1_rotated, y = X2_rotated, xend = (X1_rotated + X1_target)/2, yend = (X2_rotated + X2_target)/2), 
                 arrow = arrow(length = unit(0, 'cm')), color = "#9BBB59", linewidth = .4) +
    geom_segment(aes(x = (X1_rotated + X1_target)/2, y = (X2_rotated + X2_target)/2, xend = X1_target, yend = X2_target),
                 arrow = arrow(length = unit(0.15, 'cm')), color = "#957DB1", size = .4) +
    geom_point(aes(X1_rotated, X2_rotated), color = "#9BBB59", size = 1.6, shape = 16) + # 旋转坐标点
    geom_point(aes(X1_target, X2_target), color = "#957DB1", size = 1.6, shape = 16) + # 目的坐标点
    labs(x = 'Dim 1', y = 'Dim 2', subtitle = label, title = "Correlation analysis by Procrustes analysis") +
    geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = .4) +
    geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = .4) +
    geom_abline(intercept = 0, slope = proc_coord[1,2]/proc_coord[1,1], size = .4) +
    geom_abline(intercept = 0, slope = proc_coord[2,2]/proc_coord[2,1], size = .4) +
    theme_bw() +
    theme(axis.ticks = element_line(linewidth = .4, color = "black"),
          axis.title = element_text(size = 8, color = "black"),
          axis.text = element_text(size = 8, color = "black"),
          axis.line = element_blank(),
          plot.title = element_text(size = 10, color = "black"),
          plot.subtitle = element_text(size = 10, color = "black"),
          panel.border = element_rect(linewidth = .4, color = "black"),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          legend.text = element_text(size = 8, color = "black"),
          legend.title = element_text(size = 8, color = "black"),
          aspect.ratio = 3/4)
  message("  ggsave(file = \"procrustes_scatterplot.pdf\", width = 6, height = 4.5)\n")
  return(p)
}
