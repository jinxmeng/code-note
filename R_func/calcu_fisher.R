# Encoding: utf-8
# Author: Jinxin Meng
# Email: mengjx855@163.com
# Created Data：2023-10-28
# Modified Data: 2024-03-08
# Version: 0.1

library(tidyr)
library(dplyr)
library(tibble)

#### calcu_fisher ####
# Jinxin Meng, 20240308
# dat format:
# name  x_pos   y_pos   x_neg   y_neg
# OTU1  10      20      54      50
# OTU2  30      21      24      40
# OTU3  ..      ..      ..      ..
calcu_fisher <- function(dat, dat_colnames = NULL, plab = F, padj = F, 
                         enriched = F, enriched_by = "pval", cutoff = 0.05,
                         x = NULL, y = NULL) {
  if (!all(c("name", "x_pos", "y_pos", "x_neg", "y_neg") %in% colnames(dat)) & is.null(dat_colnames)) stop("dat field (name|x_pos|y_pos|x_neg|y_neg)")
  if (!is.null(dat_colnames)) taxonomy <- dplyr::rename(data.frame(dat, check.names = F), all_of(dat_colnames))
  if (isTRUE(enriched) & enriched_by == "padj") padj <- T
  
  dat <- column_to_rownames(dat, var = "name")
  vec_name <- rownames(dat)
  occurence <- mutate(dat, x_occur = x_pos / (x_pos + x_neg), y_occur = y_pos / (y_pos + y_neg)) %>% 
    rownames_to_column(var = "name") %>% select(name, x_occur, y_occur)
  
  pb <- txtProgressBar(style = 3, char = "#", width = 50)
  diff <- data.frame(name = character(), pval = numeric())
  for (i in 1:nrow(dat)) {
    fisher <- fisher.test(matrix(as.numeric(unlist(dat[i,])), nrow = 2, ncol = 2, byrow = T))
    diff <- add_row(diff, name = vec_name[i], pval = fisher$p.value)
    setTxtProgressBar(pb, i/nrow(dat))
  }
  close(pb)
  diff <- merge(occurence, diff, by = "name")

  if (isTRUE(plab)) diff <- add_column(diff, plab = ifelse(diff$pval < 0.001, "***", ifelse(diff$pval < 0.01, "**", ifelse(diff$pval < 0.05, "*", ""))), .after = "pval")
  if (isTRUE(padj) & isFALSE(plab)) diff <- add_column(diff, padj = p.adjust(diff$pval, method = "BH"), .after = "pval")
  if (isTRUE(padj) & isTRUE(plab)) diff <- add_column(diff, padj = p.adjust(diff$pval, method = "BH"), .after = "plab")
  if (isTRUE(enriched) & enriched_by == "pval") diff <- add_column(diff, enriched = ifelse(diff$pval < cutoff & diff$x_occur > diff$y_occur, "x", ifelse(diff$pval < cutoff & diff$x_occur < diff$y_occur, "y", "none")), .after = "pval")
  if (isTRUE(enriched) & enriched_by == "padj") diff <- add_column(diff, enriched = ifelse(diff$padj < cutoff & diff$x_occur > diff$y_occur, "x", ifelse(diff$padj < cutoff & diff$x_occur < diff$y_occur, "y", "none")), .after = "padj")
  
  if (!is.null(x)) colnames(diff)[colnames(diff)=="x_occur"] <- paste0(x, "_occur")
  if (!is.null(y)) colnames(diff)[colnames(diff)=="y_occur"] <- paste0(y, "_occur")
  if ("enriched" %in% colnames(diff)) diff <- mutate(diff, enriched = ifelse(enriched == "x", x, ifelse(enriched == "y", y, enriched)))
  
  return(diff)
}