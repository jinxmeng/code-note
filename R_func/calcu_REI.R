# Encoding: utf-8
# Source: https://github.com/twbattaglia/MicrobeDS
# Author: Johnny Hong, Ulas Karaoz, Perry de Valpine and William Fithian
# Paper: To rarefy or not to rarefy: robustness and efficiency trade-offs of rarefying microbiome data.
# DOI: https://doi.org/10.1093/bioinformatics/btac127
# Modifier: Jinxin Meng
# Created Data：2023-10-27
# Modified Data: 2023-10-27
# Version: 1.0

# REI：Rarefaction efficiency indices, if REI > 0.7, profile can be rarified.
# otu_matrix: input a otu table, taxa are column
# rarefied_depth: depth
# group: a vector that corresponding to sample
# How to explain the REI result?
# we observe that whenever the sample REI is small, the gap between the power curves is large, regardless of sample
# sizes. When the sample REI is around 0.9, the power curves essentially overlap, meaning that the effect of rarefaction on
# sensitivity is negligible; when the sample REI is around 0.7, the gap between power curves is visible but small. In practice, the default
# threshold for the sample REI can be set to 0.7; if the sample REI is below 0.7, one must beware of the drop in sensitivity due to
# rarefaction.

calu_REI <- function(otu_matrix, rarefied_depth, group, taxa_are_rows = FALSE) {
  # Compute rarefaction efficiency indices based on a given OTU matrix
  # and the rarefied depth.
  #
  # Input(s):
  #  - otu_matrix: a matrix of counts (taxa are columns)
  #  - rarefied_depth: a scalar indicating the rarefied depth
  #  - group: a vector indicating the group memberships
  #  - taxa_are_rows: a boolean, TRUE if each row represents a taxon
  #
  # Returns:
  #  - a vector of REIs
  if (taxa_are_rows) {
    otu_matrix <- t(otu_matrix)
  }
  library_size <- rowSums(otu_matrix)
  
  # Drop all observations with library sizes less than rarefied depth.
  otu_matrix <- otu_matrix[library_size >= rarefied_depth, ]
  group <- group[library_size >= rarefied_depth]
  library_size <- library_size[library_size >= rarefied_depth]
  group_lvl <- unique(group)
  n1 <- sum(group == group_lvl[1])
  n2 <- sum(group == group_lvl[2])
  
  var_prop1 <- apply(otu_matrix[group == group_lvl[1], ] / 
                       library_size[group == group_lvl[1]], 
                     2, function(col) {
                       var(col)
                     })
  var_prop_rrf1 <- apply(otu_matrix[group == group_lvl[1], ] / 
                           library_size[group == group_lvl[1]], 
                         2, function(col) {
                           var(col) + 
                             1 / rarefied_depth * mean(col * (1 - col) * (library_size[group == group_lvl[1]] - rarefied_depth) / 
                                                         (library_size[group == group_lvl[1]] - 1))
                         })
  
  var_prop2 <- apply(otu_matrix[group == group_lvl[2], ] / 
                       library_size[group == group_lvl[2]], 
                     2, function(col) {
                       var(col)
                     })
  var_prop_rrf2 <- apply(otu_matrix[group == group_lvl[2], ]  / 
                           library_size[group == group_lvl[2]], 
                         2, function(col) {
                           var(col) + 
                             1 / rarefied_depth * mean(col * (1 - col) * (library_size[group == group_lvl[2]] - rarefied_depth) / 
                                                         (library_size[group == group_lvl[2]] - 1))
                         })
  
  (var_prop1 / n1 + var_prop2 / n2) / 
    (var_prop_rrf1 / n1 + var_prop_rrf2 / n2)
}