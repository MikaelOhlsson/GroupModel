#Jaccard distance calculation
library(tidyverse)

J_m <- matrix(nrow = 26, ncol = 26, data = NA)
dimnames(J_m) <- list(subnames,subnames)

for(i in 1:26) {
  for(j in 1:26) {
    #namelist <- list()
    tmp_df <- data_frame(sp=unique(g_df$TROPHOSPECIES)) %>%
      left_join(., filter(g_df[,c(8,3,2)], WEB==subnames[i]), 
                by=c("sp" = "TROPHOSPECIES")) %>%
      left_join(., filter(g_df[,c(8,3,2)], WEB==subnames[j]), 
                suffix = c(".x", ".y"), by=c("sp" = "TROPHOSPECIES")) %>%
      #na.omit()
      filter(!is.na(Group.x) | !is.na(Group.y)) #replace with na.omit() for only comparing overlapping species
    
    glist <- unique(tmp_df$Group.x) %>% .[!is.na(.)] #Group numbers available in the "origin" web
    clist <- unique(tmp_df$Group.y) %>% .[!is.na(.)]#Comparison group numbers available in compared web
    Ja <- vector() #Vector list to save the Jaccard index for each group per network comparison, of which the mean is added to the J_m matrix for each network_x_network comparison
  
        for(u in 1:length(glist)) {
      onames <- tmp_df %>%
        filter(Group.x == glist[u]) %>%
        pull(sp)
      
      biggest <- 0
      for(v in 1:length(clist)) { #Testing all groups in the compared network
        cnames <-  tmp_df %>%
          filter(Group.y == clist[v]) %>%
          pull(sp)
        
      if(length(onames[!is.na(match(onames,cnames))]) >= biggest){ #Updates the "biggest" value IF it is indeed bigger.
        biggest <- length(onames[!is.na(match(onames,cnames))]) #Sets a new "biggest" based on intersecting overlapping species
        total <- length(unique(c(onames, cnames))) #The total number of unique species (including non-overlapping) in the two groups
        Ja[u] <- (1 - (biggest / total)) #Calculates the Jaccard index of disimmilarity for each group 
        }
      }
    }
    J_m[i,j] <- mean(Ja) #Mean of all groups put into matrix
  }
}

#Jaccard distance before normalizing
write.table(J_m, file = "../results/data/J_m.txt", sep = " ", dec = ".", quote = F)

# Symmetrical Jaccard matrix
# i.e., mean distance of A->B and B->A
 J_s <- J_m

for(i in 1:nrow(J_s)) {
  for(j in 1:ncol(J_s)) {
    if(i != j){
      J_s[i,j] <- (J_m[i,j] + J_m[j,i]) / 2
      J_s[j,i] <- (J_m[i,j] + J_m[j,i]) / 2
    }
  }
}
 
write.table(J_s, file = "../results/data/J_t.txt", sep = " ", dec = ".", quote = F)
   