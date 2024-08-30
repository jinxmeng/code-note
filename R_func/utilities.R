# info ----
# encoding: utf-8
# author: Jinxin Meng
# e-mail: mengjx855@163.com
# created data：2023-06-10
# modified data: 2024-05-20
# version: 0.1

# 2023-10-25: add function: check_file_name(), use to rename file. mainly for group file for microbiome analysis.
# 2023-11-01: update function: get_freq

library(dplyr)
library(tidyr)
library(tibble)

# 筛选一个向量中出现元素次数最多的N个元素 get_element ----
# vec: 包含所有要统计元素的向量
# 过滤标准至多只能是一个参数，不设置参数输出全部结果不进行任何过滤，过滤参数需要输入数值型的值
# 默认删除向量的NA值
# 如果想要输出一个计数表格，指定out_df为T
# 如果输出一个列表，默认输出是字符串。可选择fct和num
get_freq <- function(vec, top_n = NULL, tail_n = NULL, gt = NULL, 
                        gtn = NULL, eq = NULL, ne = NULL, lt = NULL, 
                        ltn = NULL, na.rm = T, out_df = F, out_fmt = "chr") {
  args <- c(top_n, tail_n, gt, gtn, eq, ne, lt, ltn)
  if (length(args) > 1) stop("parameter error.")
  if (isTRUE(na.rm)) vec <- vec[!is.na(vec)]
  dat <- data.frame(table(vec)) %>% arrange(desc(Freq))
  if (!is.null(top_n) & is.numeric(top_n)) { res <- head(dat, n = top_n) } 
  else if (!is.null(tail_n) & is.numeric(tail_n)) { res <- tail(dat, n = top_n) }
  else if (!is.null(gt) & is.numeric(gt)) { res <- filter(dat, Freq > gt) }
  else if (!is.null(gtn) & is.numeric(gtn)) { res <- filter(dat, Freq >= gtn) }
  else if (!is.null(eq) & is.numeric(eq)) { res <- filter(dat, Freq == eq) }
  else if (!is.null(ne) & is.numeric(ne)) { res <- filter(dat, Freq != ne) }
  else if (!is.null(lt) & is.numeric(lt)) { res <- filter(dat, Freq < lt) }
  else if (!is.null(ltn) & is.numeric(ltn)) { res <- filter(dat, Freq <= ltn) }
  else { res <- dat }
  res <- dplyr::rename(res, name = vec, freq = Freq) %>% mutate(name = as.character(name))
  if (isFALSE(out_df)) {
    res <- unlist(dplyr::select(res, name), use.names = F)
    if (out_fmt == "chr") { res <- as.character(res) }
    else if (out_fmt == "num") {res <- as.numeric(res) }
  } 
  return(res)
}

# 多位数向上或者向下取整 ----
floor_n <- function(x, n = 2) { 
  y <- x - as.numeric(str_sub(as.character(x), start = -n))
  return(y)
}
ceiling_n <- function(x, n = 2) {
  y <- x + 10^n - as.numeric((str_sub(as.character(x), start = -n)))
  return(y)
}

# 多个数据集的交集 ----
# 输入一个list，list下每个元素是一个向量
intersect_m <- function(list) {
  var <- list[[1]]
  for (i in 1:(length(list)-1)) {
    var <- intersect(var, list[[i + 1]]) 
  }
  return(var)
}

# stat_vec
tools_stat_vec <- function(x, top_n = NULL, other_name = "Other", out_vec = F) {
  dat <- as.character(x) %>% table() %>% as.data.frame(stringsAsFactors = F) %>% dplyr::rename(name = 1, n = 2)
  if (!is.null(top_n) & is.numeric(top_n)) {
    vec <- dat %>% arrange(desc(n)) %>% head(n = top_n - 1) %>% select(name) %>% unlist(use.names = F)
    dat <- dat %>% mutate(name = ifelse(name %in% vec, name, other_name)) %>% group_by(name) %>% summarise(n = sum(n)) %>% 
      arrange(desc(n))
  }
  if (isTRUE(out_vec)) dat <- dat$name
  return(dat)
}




  