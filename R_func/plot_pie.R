#### info ####
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-11-13
# modified date: 2024-03-09
# version: 0.1

library(dplyr)
library(tidyr)
library(ggplot2)

#### plot_pie ####
# dat: a data_frame contain two field, name|n;
# names: rename dat, colnames must be name and n;
# desc: descending count;
# order_list: manual list;
# add_n: add count to label;
# add_perc: add percentage to label;
# lab_cir: label circle layout;
# color: border color;
# fill; color fill for each part; if auto, fill multiple colors; if hue, fill hue palette. 
plot_pie <- function(dat, dat_colnames = NULL, desc = T, order_list = NULL, top_n = NULL, 
                     add_n = F, add_perc = T, lab_cir = F, lab_cir_flip = F, title = NULL, 
                     color = "white", fill = "auto", font_size = 2, hemi = F, start = 0) {
  if (!all(colnames(dat) %in% c("name", "n")) & is.null(dat_colnames)) stop("dat field (name|n)")
  if (!is.null(dat_colnames)) dat <- dplyr::rename(data.frame(dat, check.names = F), all_of(dat_colnames))
  
  if (!is.null(order_list)) dat <- mutate(dat, name = factor(name, order_list)) %>% arrange(name)
  if (is.null(order_list) & isTRUE(desc)) dat <- arrange(dat, desc(n))
  if (is.null(order_list) & isFALSE(desc)) dat <- arrange(dat, n)
  
  if (!is.null(top_n) & is.numeric(top_n)) {
    tmp_vec <- head(dat, top_n - 1) %>% dplyr::select(name) %>% unlist() %>% as.character()
    dat <- mutate(dat, name = as.character(name), name = ifelse(name %in% tmp_vec, name, "Other")) %>% 
      group_by(name) %>% summarise(n = sum(n)) %>% ungroup() %>% mutate(n = round(n, 2)) }
  
  if (isTRUE(desc)) dat <- arrange(dat, desc(n))
  if (isFALSE(desc)) dat <- arrange(dat, n)

  dat <- mutate(dat, perc = n/sum(n), 
                ypos = cumsum(perc) - 0.5 * perc,
                angle = ifelse(ypos < .5, 360 * ypos + 180, 360 * ypos))
  
  if (isFALSE(add_n) & isFALSE(add_perc)) dat <- mutate(dat, lab = name)
  if (isTRUE(add_n) & isTRUE(add_perc)) dat <- mutate(dat, lab = paste0(name, ", ", prettyNum(n, big.mark = ","), ", ", round(perc * 100, 1), "%"))
  if (isTRUE(add_n) & isFALSE(add_perc)) dat <- mutate(dat, lab = paste0(name, ", ", prettyNum(n, big.mark = ",")))
  if (isFALSE(add_n) & isTRUE(add_perc)) dat <- mutate(dat, lab = paste0(name, ", ", round(perc * 100, 1), "%"))
  
  if(length(fill) == 1) {
    if (fill == "auto") {
      # fill <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f", "#8dd3c7","#fdb462","#80b1d3","#fccde5","#d9ef8b","#fee391")
      fill <- c("#82ccdc","#f0ad80","#b4b3d8","#c2b75f","#87CCBA","#F9C851","#ACDEF3","#F9C6B3","#F5EAF2","#D3EDE4")
      fill <- rep(fill, time = (ceiling(nrow(dat)/10)))[1:nrow(dat)]
    } else if (fill == "hue") {
      fill <- scales::hue_pal()(nrow(dat))
    }
  } else if (length(fill) == nrow(dat)) {
    fill <- fill
  }
  
  p <- ggplot(dat, aes(x = 3, y = perc)) +
    geom_col(width = 1 , color = color, fill = fill, lwd = .4, show.legend = F) +
    coord_polar("y", start = start * pi/180) +
    theme_void() +
    theme(aspect.ratio = 1, 
          plot.title = element_text(color = "#000000", size = 8 + font_size, hjust = 0.5))
  
  if (!is.null(title)) p <- p + labs(title = as.character(title))
  if (isTRUE(lab_cir) & isFALSE(lab_cir_flip)) p <- p + geom_text(aes(x = 3.5, y = ypos, label = lab, angle = -angle), size = font_size, hjust = 0.5)
  if (isFALSE(lab_cir) & isTRUE(lab_cir_flip)) p <- p + geom_text(aes(x = 3.5, y = ypos, label = lab, angle = 270-angle), size = font_size, hjust = 0.5)
  if (isFALSE(lab_cir) & isFALSE(lab_cir_flip)) p <- p + geom_text(aes(x = 3.5, y = ypos, label = lab), size = font_size, hjust = 0.5)
  if (isTRUE(hemi)) p <- p + lims(x = c(0, 4))
  
  message("  ggsave(file = \"pie.pdf\", width = 4, height = 4)")
  
  return(p)
}