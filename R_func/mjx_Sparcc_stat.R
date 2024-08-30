####info####
# Encoding: utf-8
# Creater: 小梦
# Created date：2021-12-01
# Modified date：2022-03-08


####Function####
pacman::p_load(dplyr,tibble,vegan,tidyr)
# Function1
get.adj <- function(cor_dt, p_dt, cor = 0.6, pval = 0.05) {  # spracc输出的相关性和p值的矩阵进行合并过滤
    # 保留|cor|≥0.6 且p<0.05的值,将相关矩阵中对角线中的值（代表了自相关）转为 0
    cor_dt <- as.matrix(cor_dt)
    diag(cor_dt) <- 0
    cor_dt[abs(cor_dt) < cor] <- 0
    p_dt[p_dt < pval & p_dt >= 0] <- -1
    p_dt[p_dt >= pval] <- 0
    p_dt[p_dt==-1] <- 1
    adj <- as.matrix(cor_dt) * as.matrix(p_dt)  # 筛选后的邻接矩阵
    adj <- adj[rowSums(adj) != 0, colSums(adj) != 0]
    return(adj)
    }
  
# Function2
get.nwk.attr <- function(adj, suffix, mode = 'undirected', weighted = T, cytoscape = T, gephi = F, interaction = "all") { # 从相关性矩阵中获得网络的边和节点
    igraph <- graph_from_adjacency_matrix(as.matrix(adj), mode = 'undirected', weighted = T, diag = F)   # 邻接矩阵(adj)->igraph的邻接列表，获得含权的无向网络
    
    # 这种转换模式下，默认的边权重代表了sparcc计算的相关性（存在负值）.由于边权重通常为正值，因此最好取个绝对值，相关性重新复制一列作为记录
    E(igraph)$sparcc <- E(igraph)$weight
    E(igraph)$weight <- abs(E(igraph)$weight)
    
    if (gephi == T) { # 输出graphml格式，可使用gephi软件打开并进行可视化编辑
      write.graph(igraph, paste0(suffix,'_network.graphml'), format = 'graphml')
      }
    
    if (cytoscape == T) { # 输出gml格式，可使用 cytoscape 软件打开并进行可视化编辑
      write.graph(igraph, paste0(suffix,'_network.gml'), format = 'gml')
      }
   
    # 边列表和节点属性列表，也可以直接导入至gephi或cytoscape等网络可视化软件中进行编辑
    edge <- data.frame(as_edgelist(igraph))
    edge_list <- data.frame(source = edge[[1]], target = edge[[2]], weight = E(igraph)$weight, sparcc_r = E(igraph)$sparcc)
    node_list <- data.frame(nodes_id = V(igraph)$name, degree = degree(igraph))
    
    if(interaction == "all"){  # 所有的相互作用关系均保留
      write.table(edge_list, paste0(suffix,'_network_edge.txt'), sep = '\t', row.names = F, quote = F)
      write.table(node_list, paste0(suffix,'_network_node.txt'), sep = '\t', row.names = F, quote = F)
    } else if(interaction == "only") {  # 仅保留细菌和真菌的相互作用关系
      temp <- data.frame(node1 = substr(edge_list$source,1,1), node2 = substr(edge_list$target,1,1))
      edge_list <- edge_list[temp$node1 != temp$node2,]
      write.table(edge_list, paste0(suffix,'_network_edge.txt'), sep = '\t', row.names = F, quote = F)
      igraph <- graph_from_edgelist(as.matrix(edge_list[,1:2]), directed = F)
      node_list <- data.frame(nodes_id = V(igraph)$name, degree = degree(igraph))
      write.table(node_list, paste0(suffix,'_network_node.txt'), sep = '\t', row.names = F, quote = F)
      }

# Function3
get.nwk.stat <- function(igraph){  # network property
    # 边数量 The size of the graph (number of edges)
    num_edges <- length(E(igraph))
    # length(curve_multiple(igraph))

    # 节点数量 Order (number of vertices) of a graph
    num_vertices <- length(V(igraph))
    # length(diversity(igraph, weights = NULL, vids = V(igraph)))
    
    # 细菌和真菌的节点数量
    num_vertices_b <- length(grep("b_", as.character(V(igraph)$name)))
    num_vertices_f <- length(grep("f_", as.character(V(igraph)$name)))
    
    # 连接数(connectance) 网络中物种之间实际发生的相互作用数之和（连接数之和）占总的潜在相互作用数（连接数）的比例，可以反映网络的复杂程度
    # 同 graph.density;loops如果为TRUE,允许自身环（self loops即A--A或B--B）的存在
    connectance <- edge_density(igraph, loops = F)
    
    # 平均度(Average degree)
    # 或者为2M/N,其中M和N分别表示网络的边数和节点数
    average_degree <- mean(degree(igraph))
    
    # 正相关
    pos <- sum(E(igraph)$weight > 0)
    
    # 负相关
    neg <- sum(E(igraph)$weight < 0)
    
    # 网络统计
    nwk_stat <- data.frame(num_edge = num_edges, 
                           Positive = pos, 
                           Negative = neg,
                           num_vertice = num_vertices,
                           Bacteria = num_vertices_b,
                           Fungi = num_vertices_f,
                           Connectance = connectance, 
                           AverageDegree = average_degree,
                           row.names = )
    return(nwk_stat)
    }


####Option####
pacman::p_load(vegan,dplyr,tibble,ggplot2,tidyr,igraph)
options(stringsAsFactors = F, digits = 4)
setwd("")


####矩阵合并过滤数据####
phase <- c("p1","p2","p3")
for (i in phase) {
    cor_dt <- read.delim(paste0("Sparcc_Result/sparcc_cor_",i,".txt"), row.names = 1, sep = '\t', check.names = F)
    p_dt <- read.delim(paste0("Sparcc_Result/sparcc_pvals_",i,".txt"), row.names = 1, sep = '\t', check.names = F)
    adj <- get.adj(cor_dt, p_dt, cor = 0.6, pval = 0.05) # 过滤
    write.table(data.frame(adj, check.names = F), paste0(i, '_network_adj.txt'), sep = '\t', col.names = NA, quote = F)
  }


####网络属性文件####
phase <- c("p1","p2","p3")
for (i in phase) {  
  network_adj <- read.delim(paste0("Network_Files/", i, "_network_adj.txt"), row.names = 1, sep = '\t', check.names = F)
  get.nwk.attr(network_adj, suffix = i, cytoscape = T, interaction = "all")
  }


####网络信息统计####
phase <- c("p1","p2","p3")
nwk_stat <- rbind()
for (i in phase) {
  network_adj <- read.delim(paste0("Network_Files/", i, '_network_adj.txt'), row.names = 1, sep = '\t', check.names = FALSE)    # 网络文件
  igraph <- graph_from_adjacency_matrix(as.matrix(network_adj), mode = 'undirected', weighted = T, diag = F)    # 构建网络
  nwk_stat_i <- get.nwk.stat(igraph)    # 网络统计
  nwk_stat <- rbind(nwk_stat,nwk_stat_i) 
}

rownames(nwk_stat) <- gsub("p", "Phase_",phase)
nwk_stat <- rownames_to_column(nwk_stat, var = "phase")
write.table(nwk_stat, "Network_Files/nwk_stat.txt", sep = '\t', row.names = F, quote = F)



####可视化网络####
phase <- c("p1","p2","p3")
for (i in phase) {
  network_adj <- read.delim(paste0(i, "_network_adj.txt"), row.names = 1, sep = '\t', check.names = F)
  igraph <- graph_from_adjacency_matrix(as.matrix(network_adj), mode = 'undirected', weighted = T, diag = F)  # 构建网路 
  E_color <- E(igraph)$weight  # 修改边 
  E_color <- ifelse(E_color > 0, "#d7191c88", ifelse(E_color < 0, "#2b83ba88","grey"))  # 正相关设定为red, 负相关设定为blue 
  E(igraph)$color <- as.character(E_color)
  E(igraph)$width <- abs(E(igraph)$weight)  # 将相关系数与边的宽度关联
  E(igraph)$width_1 <- (abs(E(igraph)$weight) * 2) ^ 3  # 将相关系数与边的宽度关联
  V_color <- c()  # 修改点
  V_color[grep("b_", as.character(V(igraph)$name))] <- "#e76c44"   # 细菌设置颜色
  V_color[grep("f_", as.character(V(igraph)$name))] <- "#51bab1"   # 真菌设置颜色
  V(igraph)$color <- as.character(V_color)
  
  set.seed(2022)
  pdf(file = paste0("plot_", i, "_network.pdf"), height = 8, width = 14)
  plot(igraph, 
       main = "Co-occurrence network",
       vertex.frame.color = "black", 
       vertex.label = substr(V(igraph)$name, start = 3, stop = nchar(V(igraph)$name)), 
       vertex.label.color = "black",
       vertex.label.cex = 1.2,
       vertex.size = 18,
       layout = layout_on_sphere, 
       edge.width = E(igraph)$width_1,
       edge.lty = 1, 
       edge.curved = T, 
       margin = c(0, 0, 0, 0))
  
  legend("topright",
         bty = "n",
         legend = round(summary(E(igraph)$width)[c(1,4,6)], digits = 1),
         col = "black",
         ncol = 1, 
         cex = 1, 
         lwd = round(summary(E(igraph)$width_1)[c(1,4,6)], digits = 1), 
         text.font = 0.5, 
         text.col = "black")
  legend("right",
         bty = "n",
         legend = c("Positive", "Negative"),
         col = c("#d7191c", "#2b83ba"),
         ncol = 1, 
         cex = 1, 
         lwd = 1, 
         text.font = 0.5, 
         text.col = "black")
  legend("bottomright",
         bty = "n",
         legend = c("Bacteria", "Fungi"),
         pch = c(21,21),
         col = c("black","black"),
         pt.bg =  c("#e76c44","#51bab1"),
         pt.cex = 1.5,
         ncol = 1, 
         cex = 1,
         text.font = 0.5, 
         text.col = "black")
  dev.off()
  }


####可视化统计信息####
nwk_stat <- read.delim("Network_Files/nwk_stat.txt", sep = "\t")

# edge
edge <- nwk_stat %>% select(phase, Positive, Negative) %>% gather(key = "direction", value = "count", -phase)
ggplot(edge, aes(phase, count, fill = factor(direction, levels = c("Positive", "Negative")))) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.8) +
  geom_text(aes(label = count, y = count + 2), position = position_dodge(0.85), size = 2) +
  scale_fill_manual(values = c("#d53f5088", "#05607588")) +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "EdgeCount", fill = "Direction") + 
  theme_classic() + 
  theme(text = element_text(size = 9))
ggsave("plot_edge.pdf", width = 4, height = 3)

# node
node <- nwk_stat %>% select(phase, Bacteria, Fungi) %>% gather(key = "Class", value = "count", -phase)
ggplot(node, aes(phase, count, fill = Class)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.8) +
  geom_text(aes(label = count, y = count + 1), position = position_dodge(0.9), size = 2) +
  scale_fill_manual(values = c("#e76c44", "#51bab1")) +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "NodeCount", fill = "Class") + 
  theme_classic() + 
  theme(text = element_text(size = 9))
ggsave("plot_node.pdf", width = 4, height = 3)

# connectance
connectance <- nwk_stat %>% select(phase, Connectance)
ggplot(connectance, aes(phase, Connectance)) +
  geom_bar(stat = "identity", fill = "#DB9F75", width = 0.7) +
  geom_text(aes(label = round(Connectance, digits = 3), y = Connectance + 0.005), size = 2) + 
  scale_x_discrete(expand = c(0.2,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "Connectance") + 
  theme_classic() + 
  theme(text = element_text(size = 9))
ggsave("plot_Connectance.pdf", width = 3, height = 3)


# Average degree
average_degree <- nwk_stat %>% select(phase, AverageDegree)
ggplot(average_degree, aes(phase, AverageDegree)) +
  geom_bar(stat = "identity", fill = "#82D726", width = 0.7) +
  geom_text(aes(label = round(AverageDegree, digits = 3), y = AverageDegree + 0.15), size = 2) + 
  scale_x_discrete(expand = c(0.2,0)) + 
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "", y = "AverageDegree") + 
  theme_classic() + 
  theme(text = element_text(size = 9))
ggsave("plot_AverageDegree.pdf", width = 3, height = 3)

