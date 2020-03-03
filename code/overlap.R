#Species overlap between the different subregions as well as the metaweb
library(tidyverse)
library(igraph)

#Extract species lists as ordered in adjacency matrices for later use when linking species to group numbers
sp_net <- list()
for(n in 1:26) {
  x <- as.matrix(get.adjacency(net[[n]]))
  sp_net[[n]] <- colnames(x) #List of species lists for each subnet
}

#Prepare matrix
overlap_matrix <- matrix(nrow=26, ncol=26, data = rep(0, 26*26))
dimnames(overlap_matrix) <- list(subnames,subnames)

#...and do the comparison                         
for(i in 1:length(sp_net)) {
  for(j in 1:length(sp_net)) {
  overlap_matrix[i,j] <- sum(!is.na(match(sp_net[[i]], sp_net[[j]]))) / length(unique(c(sp_net[[i]], sp_net[[j]]))) 
  }
}

write.table(overlap_matrix, file = "../results/data/overlap_matrix.txt", sep = " ", dec = ".", quote = F)