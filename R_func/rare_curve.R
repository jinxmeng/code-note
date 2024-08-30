# info ----
# Encoding: utf-8
# Author: Jinxin Meng
# Email: mengjx855@163.com
# Created Data：2022-05-29
# Modified Data: 2023-12-10
# Version: 1.1

library(dplyr)
library(tibble)
library(tidyr)
library(vegan)
library(ggplot2)

# calcu_specaccum ----
# Rarefaction curve analysis | accumulative species
# calcu_specaccum <- function(otu, step = 5, permutations = 99) {
#   dat <- t(otu)
#   n <- dim(dat)[1]
#   seq <- seq(1, n, step)
#   if (max(seq) != n) seq <- c(seq, n)
#   result <- data.frame(name = numeric(), obs = numeric(), sd = numeric())
#   for (i in seq) {
#     vec <- replicate(permutations, ifelse(i == 1,
#                                           sum(dat[sample(n, i, replace = F),] > 0),
#                                           sum(colSums(dat[sample(n, i, replace = F),]) > 0)))
#     result <- add_row(result, name = i, obs = mean(vec), sd = sd(vec))
#   }
#   return(result)
# }

# calcu_specaccum ----
# Rarefaction curve analysis | accumulative species
calcu_specaccum <- function(otu, permutations = 10, method = "random"){
  sampling <- specaccum(t(otu), method = method, permutations = permutations)
  result <- data.frame(name = sampling$sites, obs = sampling$richness, sd = sampling$sd)
  return(result)
}

#### calcu_specaccum_by_group ####
# Rarefaction curve analysis | accumulative species by group
calcu_specaccum_by_group <- function(otu, group, group_colnames = NULL, group_order = NULL,
                                     permutations = 10, method = "random"){
  if (!all(c("sample", "group") %in% colnames(group)) & is.null(group_colnames)) stop("group field (sample|group)")
  if (!is.null(group_colnames)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_colnames))
  group <- filter(group, sample %in% colnames(otu))
  if (is.null(group_order)) group_order <- unique(group$group)
  otu <- t(otu)
  result <- rbind()
  for (i in group_order) {
    sample_i <- subset(group, group %in% i) %>% select(sample) %>% unlist(use.names = F)
    otu_i <- otu %>% subset(rownames(.) %in% sample_i)
    sampling <- specaccum(otu_i, method = method, permutations = permutations)
    sampling <- data.frame(name = sampling$sites, obs = sampling$richness, sd = sampling$sd) %>% add_column(group = i)
    result <- rbind(result, sampling)
  }
  return(result)
}

#### calcu_specaccum_by_depth ####
# Rarefaction curve analysis | depth
# calcu_specaccum_by_depth <- function(otu, step = 10000000) {
#   otu <- t(otu)
#   n <- dim(otu)[1]
#   dep <- rowSums(otu)
#   # 按照测序量抽平数据
#   sampling <- seq(0, max(dep), step)
#   if (max(sampling) < max(dep)) sampling <- c(sampling, max(dep))
#   dat <- rbind()
#   for (i in sampling) {
#     rarefied_otu <- rrarefy(otu, i)
#     dat_i <- estimateR(rarefied_otu)[1,] %>% data.frame(obs = .) %>% 
#       rownames_to_column(var = "sample") %>% add_column(dep = i)
#     dat <- rbind(dat, dat_i)
#   }
#   return(dat)
# }

#### plot_specaccum ####
# dat: a data.frame with field name|obs|sd
# dat_colnames: a vector to rename dat
plot_specaccum <- function(dat, dat_colnames = NULL, color = "#de2726", add_errorbar = T, 
                           add_ribbon = F, fill = NULL, lty = "solid", aspect_ratio = 1,
                           xlab = "Cumulative samples", ylab = "Cumulative features", 
                           title = "Rarefaction curve analysis") {
  if (!all(colnames(dat) %in% c("name", "obs", "sd")) & is.null(dat_colnames)) stop("dat field (name|obs|sd)")
  if (!is.null(dat_colnames)) dat <- data.frame(dat, check.names = F) %>% dplyr::rename(all_of(dat_colnames))
  if (isTRUE(add_errorbar) & isTRUE(add_ribbon)) add_errorbar = FALSE
  if (is.null(fill)) fill = "#fcbba1"
  
  plotdat <- data.frame(dat, check.names = F)
  p <- ggplot(plotdat, aes(x = name, y = obs))
  if (isTRUE(add_errorbar)) p <- p + geom_errorbar(aes(ymin = obs - sd, ymax = obs + sd), width = .6, lwd = .4, color = color, show.legend = F)
  if (isTRUE(add_ribbon)) p <- p + geom_ribbon(aes(ymin = obs - sd, ymax = obs + sd), fill = fill, show.legend = F)
  p <- p + geom_line(lwd = .4, color = color, linetype = lty, show.legend = F) +
    labs(x = xlab, y = ylab, title = title) +
    scale_x_continuous(expand = c(.02, .02)) +
    scale_y_continuous(expand = c(.02, .02)) +
    theme_bw() +
    theme(axis.line = element_line(linewidth = .4, color = "#000000"),
          axis.ticks = element_line(linewidth = .4, color = "#000000"),
          axis.text = element_text(size = 8, color = "#000000"),
          axis.title = element_text(size = 8, color = "#000000"),
          plot.title = element_text(size = 10, color = "#000000"),
          panel.grid = element_blank(),
          aspect.ratio = aspect_ratio)
  message("  ggsave(file = \"plot.rarefaction_specaccum.pdf\", width = 4, height = 4)")
  return(p)
}


#### plot_specaccum_by_group ####
# dat: a data.frame with field name|obs|sd
# dat_colnames: a vector to rename dat
plot_specaccum_by_group <- function(dat, dat_colnames = NULL, add_errorbar = T, add_lab_to_plot = F,
                                    add_ribbon = F, fill = NULL, lty = "solid", aspect_ratio = 1,
                                    xlab = "Number of samples", ylab = "Number of features", 
                                    title = "Rarefaction curve analysis") {
  if (!all(colnames(dat) %in% c("name", "obs", "sd", "group")) & is.null(dat_colnames)) stop("dat field (name|obs|sd|group)")
  if (!is.null(dat_colnames)) dat <- data.frame(dat, check.names = F) %>% dplyr::rename(all_of(dat_colnames))
  if (isTRUE(add_errorbar) & isTRUE(add_ribbon)) add_errorbar = FALSE
  if (is.null(fill)) fill = "grey85"
  
  plotdat <- data.frame(dat, check.names = F)
  p <- ggplot(plotdat, aes(x = name, y = obs, color = group, group = group))
  if (isTRUE(add_errorbar)) p <- p + geom_errorbar(aes(ymin = obs - sd, ymax = obs + sd), width = .3, lwd = .4, show.legend = F)
  if (isTRUE(add_ribbon)) p <- p + geom_ribbon(aes(ymin = obs - sd, ymax = obs + sd), fill = fill, color = NA, show.legend = F)
  p <- p + geom_line(lwd = .4, linetype = lty) +
    labs(x = xlab, y = ylab, title = title) +
    scale_x_continuous(expand = c(.02, .02)) +
    scale_y_continuous(expand = c(.02, .02)) +
    theme_bw() +
    theme(axis.line = element_line(linewidth = .4, color = "#000000"),
          axis.ticks = element_line(linewidth = .4, color = "#000000"),
          axis.text = element_text(size = 8, color = "#000000"),
          axis.title = element_text(size = 8, color = "#000000"),
          plot.title = element_text(size = 10, color = "#000000"),
          panel.grid = element_blank(),
          aspect.ratio = aspect_ratio)
  if (isTRUE(add_lab_to_plot)) {
    plotlab <- group_by(plotdat, group) %>% group_modify(~.x[which.max(.x$name),])
    p <- p + geom_label(data = plotlab, aes(x = name, y = obs, label = group), inherit.aes = F, size = 1.2) +
      guides(color = "none")
  }
  message("  ggsave(file = \"plot.rarefaction_specaccum_by_group.pdf\", width = 4, height = 4)")
  return(p)
}
