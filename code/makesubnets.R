library(igraph)
library(readr)
library(tidyverse)

#Species list with trophic levels
splist <- read_delim("../data/kortsch/SpeciesList_Kortsch.txt", delim="\t", col_names = T)
splist <- data_frame(ABBREVIATION=splist$ABBREVIATION, TROPHOSPECIES=splist$TROPHOSPECIES, GROUP=splist$GROUP) #Reorder columns to have abbreviations first, matching the edge list later

#Subregion network df
subr <- read_delim("../data/kortsch/SubregionSpecies.txt", delim = "\t", col_names = F)
subr[1,202] <- "SEB_SPP" #Typo fix
write.table(subr, "../data/kortsch/SubregionSpecies.txt", row.names = F, col.names=F, sep = "\t", quote = F)


#Create new edge lists for each network...
net <- list()
for(n in 1:nrow(subr)){
  pwl <- read_delim("../data/kortsch/PairwiseList_seb.txt", delim = " ", col_names = T) #PairwiseList_Kortsch.txt with added Sebastes spp. links

# # # Filtering the metaweb based on subregion species compositions to generate subregion networks # # #
#Removes rows in the edge list where the *predator* is not present in the current subregion
if(n > 1){
  for(i in nrow(pwl):1) {
    if(!is.na(match(pwl[i,1], subr[1,])))  pos <- match(pwl[i,1], subr[1,])
    if(subr[n,pos] == 0)  pwl <- pwl[-i,]
  }
  
  #Removes rows in the edge list where the *prey* is not present in the current subregion
  for(i in nrow(pwl):1) {
    if(!is.na(match(pwl[i,2], subr[1,])))  pos <- match(pwl[i,2], subr[1,])
    if(subr[n,pos] == 0)  pwl <- pwl[-i,]
  }
}  
print(n-1)
#And merge the edge lists and species data into a new graph object, saved in list "net"
net[[n]] <- graph_from_data_frame(d=pwl, vertices = splist, directed=T) #Also re-adds uninteracting vertices
net[[n]] <- igraph::delete.vertices(net[[n]], igraph::degree(net[[n]])==0) #Removes non-interacting vertices
}
saveRDS(net, file="../data/graphs_seb.rds")
#DONE

#Writing subregion networks into individual files (as adjacency matrices)
#Adjacency matrices with names
for(i in 1:length(net)){
  x <- as.matrix(get.adjacency(net[[i]]))
  write.table(x, file=paste0("../data/sebwebs_with_names/", ifelse(i==1, "fullnet", subnames[i]),".txt"), row.names = T, col.names = NA, quote = F)
}

#Adjacency matrices without names (required for the group model code)
for(i in 1:length(net)){
  x <- as.matrix(get.adjacency(net[[i]]))
  write.table(x, file=paste0("../data/sebwebs/", ifelse(i==1, "fullnet", subnames[i]),".txt"), row.names = F, col.names = F, quote = F)
}
