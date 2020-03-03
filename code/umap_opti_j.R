# UMAP clustering optimization for smallest variation of the clusters
# This version optimizes the overall mean of the clusters, 
# allowing e.g. one cluster with high variability if the rest of the clusters are very low

library(umap)
library(tidyverse)
library(stringr)
library(dbscan)
library(reshape2)
library(doMC)

#Note: Also requires the python package "umap-learn"

#Comparison matrix to be tested (Jaccard distance in this case)
c_m <- as.matrix(read.delim("../results/data/J_t.txt", sep= " ", dec=".", header=T)) %>%
  .[-1,-1] #Exclude metaweb

umap.mod <- umap.defaults #Dummy umap settings to be modified
umap.mod$min_dist <- 0.001 # Low values enhances local similarities in the projection, which helps HDBSCAN

ClusterOptimization <- function(Iterations, MinCluster, Unclustered) {
  netnames <- dimnames(c_m)[[1]]
  topC <- data.frame(neighbors=NA, mindist=NA, Cmean=NA, Cclust=NA, CX=NA, CY=NA)
  Objects <- nrow(c_m) #Number of networks
  
  #Parallelizing forloops 
  top_c <- foreach(m=1:Iterations, .combine= rbind) %dopar% { #Each iteration tests all combination of UMAP scale and HDBSCAN cluster min size
    set.seed(floor(runif(1,101,1000000))) 
    allclust <- data.frame(neighbors=NA, minpts=NA, Cmean=NA, Cclust=NA, CX=NA, CY=NA) 
    count <- 0
    for(i in 2:Objects){ # i = UMAP scale, small numbers focus on more local scales
        count <- count + 1
        
        umap.mod$n_neighbors <- i #UMAP scale parameter, max is the number of networks
        
        umap_c <- umap(c_m, config = umap.mod, method = "umap-learn") # method requires install the python package "umap-learn"
        hdb_c <- hdbscan(umap_c$layout, minPts=(MinCluster))
        
       # MinDist <- MinDist + 0.001 #increasing minimum distance during the j-loop
        
        c_clusters <- data.frame(WEB=netnames, Cluster=hdb_c$cluster)

        # Creates data frames with the network pairs for each cluster
        # value column is the dissimilarity score from Jaccard index
        c_melt <- left_join(c_clusters[,c(1,2)], melt(c_m), by=c("WEB" = "Var1")) %>%
          left_join(c_clusters[,c(1,2)], by=c("Var2" = "WEB")) %>%
          filter(Cluster.x == Cluster.y & Cluster.x > 0)
        
        allclust[count,1] <- i #Setting for umap projection (projection scale)
        allclust[count,2] <- umap.mod$min_dist #UMAP minimum distance
        
        # Calculates the mean cluster dissimilarity
        allclust[count,3] <- c_melt[c_melt$value > 0,] %>% #Excludes comparing networks with themselves
          filter(Cluster.x > 0) %>% #Excludes unclustered networks
          group_by(Cluster.x) %>%
          summarise(avg=mean(value)) %>% #Mean for each cluster
          summarize(totmean=mean(avg)) #Overall mean
        
        
        #"Cheating" with collapsing cluster lists and layouts (umap coordinates) instead of nested data frames
        allclust[count,4] <- paste(hdb_c$cluster, collapse = " ", sep = " ")
        allclust[count,5] <- paste(umap_c$layout[,1], collapse = " ", sep = " ")
        allclust[count,6] <- paste(umap_c$layout[,2], collapse = " ", sep = " ")
    #  }
    }
    
    
    cfilt <- filter(allclust, str_count(allclust$Cclust, "0") <= Unclustered) #Limit the number of unclustered regions to Unclustered variable (0 for now)
    
    cfilt[match(min(cfilt[,3], na.rm = T), cfilt[,3]), c(1,2,3,4,5,6)] #Extract the best iteration based on JACCARD INDEX
    #cfilt[match(max(cfilt[,3], na.rm = T), cfilt[,3]), c(1,2,3,4,5,6)] #Extract the best iteration based on SPECIES OVERLAP
    }
  
  
  Cn <- match(min(top_c[,3]), top_c[,3]) #Best JACCARD clustering (lowest cluster average)
  #Cn <- match(max(top_c[,3]), top_c[,3]) #Best SPECIES OVERLAP clustering (highest cluster average)
  Ci <- top_c[Cn,1] # UMAP Neighbor number 
  Cj <- top_c[Cn,2] # UMAP Minimum distance
  c_score <- top_c[Cn,3]
  c_clusters <- data.frame(WEB=netnames, Cluster=strsplit(top_c[Cn,4], split = " "), 
                           X=strsplit(top_c[Cn,5], split = " "), Y=strsplit(top_c[Cn,6], split = " ")) 
  colnames(c_clusters)[1:4] <- c("WEB", "Cluster", "X", "Y")
  
  topclusters <- list(c_clusters, c_score, Ci, Cj) #Result list
  names(topclusters) <- c("c_clusters", "c_score", "c_neigh", "c_mindist")
  
  return(topclusters)
}

# # # Optimization settings # # #
registerDoMC(10) #Number of cores to use
it <- 1000 	#Iterations, 1000+ recommended (10 000 for the study)
clu <- 3 	#Minimum cluster size 
uncl <- 0 	#Number of allowed unclustered networks

# roughly 50 iterations per core and minute
jd_topclusters <- ClusterOptimization(Iterations=it,
                                   MinCluster=clu,
                                   Unclustered=uncl) 


#saveRDS(jd_topclusters, paste0("../results/data/clustering/clusteropti_",it,"i_",clu,"c_", uncl,"u_jd.rds"))
