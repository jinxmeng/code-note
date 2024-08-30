#### info ####
# encoding: utf-8
# author: Jinxin Meng
# e-mail: mengjx855@163.com
# created data：2024-03-04
# modified data: 2024-03-04

library(paletteer)

#### colors_pick ####
# colors_pick
# name: 预设的颜色集
# n: 输出颜色的数量
# 不指定参数默认输出颜色集的名称
colors_pick = function(name = NULL, n = NULL) {
  colors_set = c("set2_8","set3_12","accent_8","paired_12","paired2_12","paired_20","hue_20","comparison_3", 
             "comparison2_2","comparison3_2","continuous_ryb_7","continuous_pwg_7","continuous_rwb_7","continuous_hwc_7")
  
  if (is.null(name) & is.null(n)) return(paste0("Option: ", paste(colors_set, collapse = ", ")))
  
  colors_list <- list(
    set2_8 = c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494","#b3b3b3"),
    set3_12 = c("#80b1d3","#b3de69","#fdb462","#8dd3c7","#bc80bd","#fb8072","#ffed6f","#fccde5","#bebada","#ccebc5","#ffffb3","#d9d9d9"),
    accent_8 = c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17","#666666"),
    paired_12 = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928"),
    paired2_12 = c("#1f78b4","#a6cee3","#33a02c","#b2df8a","#e31a1c","#fb9a99","#ff7f00","#fdbf6f","#6a3d9a","#cab2d6","#b15928","#ffff99"),
    paired_20 = paletteer::paletteer_d("ggthemes::Classic_20"),
    hue_20 = c("#f8766d","#ea8331","#d89000","#c09b00","#a3a500","#7cae00","#39b600","#00bb4e","#00bf7d","#00c1a3","#00bfc4","#00bae0","#00b0f6","#35a2ff","#9590ff","#c77cff","#e76bf3","#fa62db","#ff62bc","#ff6a98"),
    comparison_3 = c("#56B4E9", "#E69F00", "#999999"),
    comparison2_2 = c("#f46d43","#74add1"),  # 多种免疫病病毒组研究中病例和对照组用的颜色
    comparison3_2 = c("#f77d4d", "#1f78b4"),  # MDD,IBS,ACVD病例对照中的颜色方案
    continuous_ryb_7 = c("#d53e4f","#fc8d59","#fee08b","#ffffbf","#e6f598","#99d594","#3288bd"), # 渐变色 红-黄-蓝
    continuous_pwg_7 = c("#c51b7d","#e9a3c9","#fde0ef","#f7f7f7","#e6f5d0","#a1d76a","#4d9221"), # 渐变色 紫-白-绿
    continuous_rwb_7 = c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac"), # 渐变色 红-白-蓝
    continuous_hwc_7 = c("#8c510a","#d8b365","#f6e8c3","#f5f5f5","#c7eae5","#5ab4ac","#01665e") # 渐变色 褐-白-青
  )
  
  if (!is_null(name) & !name %in% colors_set) stop("colors set name error. please using color_pick() to check")
  
  if (!is.null(name) & is.null(n)) {
    x = colors_list[[name]]
  } else {
    x = colors_list[[name]]
    len = length(x)
    x <- rep(x, times = ceiling(n/len))[1:n]
  }
  return(x)
}

#### colors_continuous ####
# colors_continuous
# name: 预设的颜色集
# n: 输出颜色的数量
# 不指定参数默认输出颜色集的名称
colors_continuous = function(name = NULL, n = NULL) {
  colors_set = c("hue","yh","bwr","hwc","ywp","rwb","ywb","rwg")
  if (is.null(name) & is.null(n)) return(paste0("Option: ", paste(colors_set, collapse = ", ")))
  if (is.null(name) & !is.null(n)) stop("Both name and n parameters must be specified simultaneously")
  if (!is.null(name) & is.null(n)) stop("Both name and n parameters must be specified simultaneously")
  
  if (name == "hue") {
    x = scales::hue_pal()(n)
  } else if (name == "yh" ) {
    x = colorRampPalette(c("#F3DA7D","#A13043"))(n) # 渐变色 黄-褐
  } else if (name == "bwr" ) {
    x = colorRampPalette(c("#3288bd", "#ffffff", "#d53e4f"))(n) # 渐变色 蓝-白-红
  } else if (name == "hwc" ) {
    x = colorRampPalette(c("#a6611a","#dfc27d","#f5f5f5","#80cdc1","#018571"))(n) # 渐变色 褐-白-青
  } else if (name == "ywp" ) {
    x = colorRampPalette(c("#e66101","#fdb863","#f7f7f7","#b2abd2","#5e3c99"))(n) # 渐变色 黄-白-紫
  } else if (name == "rwb" ) {
    x = colorRampPalette(c("#ca0020","#f4a582","#f7f7f7","#92c5de","#0571b0"))(n) # 渐变色 红-白-蓝
  } else if (name == "ywb" ) {
    x = colorRampPalette(c("#d7191c","#fdae61","#ffffbf","#abd9e9","#2c7bb6"))(n) # 渐变色 红-黄-蓝
  } else if (name == "rwg" ) {
    x = colorRampPalette(c("#f46d43","#fee08b","#ffffff","#d9ef8b","#66bd63"))(n) # 渐变色 红-白-绿
  } 
  return(x)
}


