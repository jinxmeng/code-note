# Jinxin Meng, 20240305, 20240712 ------------------
pacman::p_load(dplyr, tidyr, tibble, ggClusterNet, phyloseq, igraph)

# calcu_MEN ---------------
# 微生物生态网络
# 输入profile表，物种表
# metadata 对OTU的注释信息，行名和profile行名一致，随便伪造一个就行，界门纲目科属种
#       Phylum  Family ...
# OTU_1 phy_1   fam_1 ...
# OTU_2 phy_2   fam_2 ...
# OTU_3 phy_3   fam_3 ...
# ...   ...     ...   ...

calcu_MEN <- function(profile, metadata, p.threshold = 0.05, r.threshold = 0.5,
                      scale = "TMM", method = "spearman", p.adj = "BH", 
                      title = "", n.hub = F){
  otu_table <- otu_table(profile, taxa_are_rows = T)
  tax_table <- metadata %>% 
    rename_all(\(x) stringr::str_to_title(x)) %>% 
    as.matrix() %>% 
    tax_table()
  ps <- phyloseq(otu_table, tax_table)
  out <- list()
  
  # 相关性分析
  cor_obj <- cor_Big_micro2(ps, N = 0, p.threshold = p.threshold, r.threshold = r.threshold,
                            scale = scale, method = method, p.adj = p.adj) 
  r_data <- cor_obj[[1]]
  
  # 一个list，第一个元素是边数据，第二个是节点数据
  edge_obj <- nodeEdge(corr = r_data)
  igraph_obj <- igraph::graph_from_data_frame(edge_obj[[1]], directed = F, vertices = edge_obj[[2]])
  
  # 网络属性
  # net_properties (14) 
  # net_properties.2
  # net_properties.3
  # net_properties.4
  attr_data <- net_properties.2(igraph_obj, n.hub = n.hub) %>% 
    data.frame() %>% 
    rownames_to_column("name")
  
  # model_igraph2布局
  graph <- model_igraph2(cor = r_data, method = "cluster_fast_greedy", seed = 2024)
  
  # node
  otu_table <- t(vegan_otu(ps)) %>% 
    data.frame(check.names = F)
  tax_table <- vegan_tax(ps) %>% 
    data.frame(check.names = F)
  node_data <- graph[[1]] %>% 
    nodeadd(otu_table = otu_table, tax_table = tax_table)
  
  # model & color
  model <- graph[[2]]
  color_node <- data.frame(model = model$model, color = model$color) %>% 
    distinct(model, .keep_all = T) %>% 
    pull(2, name = 1)
  
  # edges
  edge_data <- edgeBuild(cor = r_data, node = graph[[1]])
  edges <- inner_join(
    select(model, OTU, model, color) %>% 
      right_join(edge_data, by = c("OTU" = "OTU_1")) %>% 
      rename(OTU_1 = OTU, model1 = model, color1 = color),
    select(model, OTU, model, color) %>%
      right_join(edge_data, by = c("OTU" = "OTU_2")) %>%
      rename(OTU_2 = OTU, model2 = model, color2 = color)) %>% 
    mutate(edge_class = ifelse(model1 == model2, model1, "across"), # 来自不同模块的连接，填充为灰色
           edge_color = ifelse(model1 == model2, color1, "#C1C1C1"))
  
  color_edge <- distinct(edges, edge_class, .keep_all = T) %>% 
    select(edge_class, edge_color) %>% 
    pull(2, name = 1)
  
  p <- ggplot() + 
    geom_segment(data = edges, aes(x = X1, y = Y1, xend = X2, yend = Y2, color = edge_class), 
                 linewidth = 1, show.legend = F) + 
    scale_colour_manual(values = color_edge) +
    ggnewscale::new_scale_color() +
    geom_point(data = model, aes(X1, X2, color = model), size = 4, show.legend = F) +
    scale_colour_manual(values = color_node) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    labs(title = title) +
    theme(axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          plot.title = element_text(size = 10, hjust = 0.5),
          legend.background = element_rect(colour = NA),
          panel.background = element_rect(fill = "white",  colour = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 1, plot.margin = unit(c(1, 1, 1, 1), "cm"))
  
  out[["profile"]] <- profile
  out[["metadata"]] <- metadata
  out[["cor_obj"]] <- cor_obj
  out[["igraph_obj"]] <- igraph
  out[["network_attr"]] <- attr_data
  out[["model_obj"]] <- model
  out[["node_data"]] <- node_data
  out[["edge_data"]] <- edge_data
  out[["gplot"]] <- p 
  return(out)
}


