library(tidyverse)
library(igraph)
library(NetIndices)
library(cheddar)
library(stringr)

net <- readRDS("../data/graphs_seb.rds") #File containing graphs for all networks (from makesubnets.R)
subr <- read_delim("../data/kortsch/SubregionSpecies.txt", delim = "\t", col_names = F)

#Extract species lists as ordered in adjacency matrices for later use when linking species to group numbers
sp_net <- list()

tl_df <- data_frame(Species=NA, TL=NA, OI=NA, Network=NA)
tl_df <- tl_df[-1,]
for(n in 1:26) {
  #Trophic level
  x <- net[[n]] %>%
    get.adjacency() %>%
    as.matrix() %>%
    TrophInd(Tij=., Dead = 1) %>%  #Tij instead of Flow to specify the direction of the matrix, from column to row, Dead = 1 specifies detritus (species 1)
    cbind(., data.frame(Network=ifelse(n==1, "Metaweb", paste0("sub", subr[n,1])))) %>%
    rownames_to_column(var="Species")
  
  #Generality
  x <- net[[n]] %>%
    get.adjacency() %>%
    as.matrix() %>%
    rowSums() %>%
    data.frame(Gen=.) %>%
    rownames_to_column(., var="Species") %>%
    left_join(x,.)
  #Vulnerability
  x <- net[[n]] %>%
    get.adjacency() %>%
    as.matrix() %>%
    colSums() %>%
    data.frame(Vul=.) %>%
    rownames_to_column(., var="Species") %>%
    left_join(x,.)
  #Total links
  x <- mutate(x, Links=Gen+Vul)
  
  tl_df <- bind_rows(tl_df, x)
}

# # # Adding groups from SGM # # #
setwd("../results/groups")
results <- list.files(path = ".", pattern = "Marginal") #Result files
webs <- list.files(path= ".", pattern = "*.txt$") #Adjacency matrix files, used for names
g_df <- data_frame(Species=NA, WEB=NA, Group=NA, TL=NA, Gen=NA, Vul=NA, Links=NA) #Data frame which all subnet groupings will be added to
g_df <- g_df[-1,]

# Metaweb grouping added first, separately
lowest <- str_sort(results[str_detect(pattern="full", string = results) == TRUE], decreasing = FALSE)[1] #Picks the grouping with best marginal
  x <- read.table(lowest, sep=" ", col.names=F, as.is = T) #Loads best subnet grouping result
  
  g <- data_frame(Species=V(net[[1]])$name, 
                  WEB="Metaweb", 
                  Group=x$FALSE.+1, 
                  TL=tl_df$TL[tl_df$Network=="Metaweb"], 
                  Gen=tl_df$Gen[tl_df$Network=="Metaweb"], 
                  Vul=tl_df$Vul[tl_df$Network=="Metaweb"], 
                  Links=tl_df$Links[tl_df$Network=="Metaweb"]) #Attach groupings to the correct species
  g_df <- bind_rows(g_df, g) #Joins the subregion to the full set while matching species names


# Subnet groupings 
# File lists (webs and results) are not ordered as the original sub region list (subr[,1]), so it takes some extra matching of strings first
for(i in 2:nrow(subr)) {
  print(i)
  lowest <- NA
  w <- str_which(pattern=as.character(paste0("^",subr[i,1])), string=str_replace_all(webs, "sub", "")) #Which webs object corresponds to the current subregion
  lowest <- str_sort(results[str_detect(pattern=webs[[w]], string = results) == TRUE], decreasing = FALSE)[1] #Picks the grouping with best marginal
  
    x <- read.table(lowest, sep=" ", col.names=F, as.is = T) #Loads best subnet grouping result
    
    g <- data_frame(Species=V(net[[i]])$name, 
                    WEB=substr(webs[[w]], 1, nchar(webs[[w]])-4), 
                    Group=x$FALSE.+1,
                    TL=tl_df$TL[tl_df$Network==substr(webs[[w]], 1, nchar(webs[[w]])-4)], 
                    Gen=tl_df$Gen[tl_df$Network==substr(webs[[w]], 1, nchar(webs[[w]])-4)], 
                    Vul=tl_df$Vul[tl_df$Network==substr(webs[[w]], 1, nchar(webs[[w]])-4)], 
                    Links=tl_df$Links[tl_df$Network==substr(webs[[w]], 1, nchar(webs[[w]])-4)]) #Attach groupings to the correct species
    g_df <- bind_rows(g_df, g) #Joins the subregion to the full set while matching species names
}
setwd("../../code")
#Adding taxonomy
tax <- read_delim("../data/kortsch/SpeciesList.txt", delim="\t")
colnames(g_df)[1] <- "ABBREVIATION"

g_df <- left_join(g_df, tax[,1:6], by="ABBREVIATION")

#Saving results
write.table(g_df, file="../results/data/g_df.txt", quote = F, row.names = F, sep = "\t")
#Done