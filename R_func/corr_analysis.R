# Encoding: utf-8
# Author: Jinxin Meng
# Email: mengjx855@163.com
# Created Data：2022-04-25
# Modified Data: 2023-10-28
# Version: 1.0

#### Function1 ####
get.adj <- function(cor_dt, p_dt, r = 0.6, p = 0.01) {  # corr.test 输出的相关性和p值的矩阵进行合并过滤
  # 保留|cor|≥0.6 且p<0.05的值,将相关矩阵中对角线中的值（代表了自相关）转为 0  
  cor_dt <- as.matrix(cor_dt)
  cor_dt[abs(cor_dt) < r] <- 0
  p_dt[p_dt >= p] <- -1
  p_dt[p_dt < p & p_dt >= 0] <- 1
  p_dt[p_dt==-1] <- 0
  adj <- as.matrix(cor_dt) * as.matrix(p_dt)  # 筛选后的邻接矩阵
  adj <- adj[rowSums(adj) != 0, colSums(adj) != 0]
  return(adj)
}

#### Function2 ####
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
}

#### Function3 ####
get.nwk.stat <- function(igraph){  # network property
  # 边数量 The size of the graph (number of edges)
  num_edges <- length(E(igraph))
  # length(curve_multiple(igraph))
  
  # 节点数量 Order (number of vertices) of a graph
  num_vertices <- length(V(igraph))
  # length(diversity(igraph, weights = NULL, vids = V(igraph)))
  
  # 细菌和真菌的节点数量
  num_vertices_b <- length(grep("^b_", as.character(V(igraph)$name)))
  num_vertices_f <- length(grep("^f_", as.character(V(igraph)$name)))
  
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
