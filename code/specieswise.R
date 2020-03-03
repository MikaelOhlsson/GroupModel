#Species-wise group turnover, calculated as the mean turnover from network A to all other networks
library(tidyverse)

#Data frame for collecting the species- and web-specific group turnover
Ti <- tibble(TROPHOSPECIES = "S", WEB = "S", Turnover = 0.0) %>% .[-1,] 

for(u in 1:length(subnames)){
  A <- g_df %>% #First web to use as base comparison
    filter(WEB == subnames[u]) %>%
    dplyr::select(TROPHOSPECIES, Group)
  TA <- tibble(TROPHOSPECIES=unique(g_df$TROPHOSPECIES))

  for(v in 1:length(subnames)){
    if(u != v){
      B <- g_df %>% #all other networks
        filter(WEB == subnames[v]) %>%
        dplyr::select(TROPHOSPECIES, Group)
      
      SharedSp <- filter(intersect(A[,1], B[,1])) %>% #only looking at shared species between A and B
        left_join(A) %>%
        left_join(B, by=c("TROPHOSPECIES"), suffix = c(".A", ".B"))  
      
      Sab <- nrow(SharedSp) #number of shared species
      
      #Generating relationship matrices for all species pairs in both networks
      Qa <- matrix(nrow = Sab, ncol = Sab, data = NA) 
      dimnames(Qa) <- list(SharedSp$TROPHOSPECIES, SharedSp$TROPHOSPECIES)
      Qb <- Qa
      
      for(i in 1:Sab) {
        for(j in 1:Sab){
          Qa[i,j] <- ifelse(SharedSp$Group.A[i] == SharedSp$Group.A[j], 1, 0)
          Qb[i,j] <- ifelse(SharedSp$Group.B[i] == SharedSp$Group.B[j], 1, 0)
        }
      }
      
      #Calculating the mean turnover for each species
      TAB <- rowSums(Qa * (1 - Qb) + Qb * (1 - Qa)) / (Sab-1)
      TAB <- tibble(TROPHOSPECIES = names(TAB), T=TAB)
      TA <- left_join(TA, TAB, by=c("TROPHOSPECIES")) #Adding columns for each web comparison to A...
    }
  }
 TA <- tibble(TROPHOSPECIES = TA$TROPHOSPECIES, WEB = subnames[u], Turnover = rowMeans(TA[,-1], na.rm = T)) #...from which the mean is derived
 
 Ti <- rbind(Ti, TA)
}

#Binding the complete results
g_df <- left_join(g_df, Ti, by=c("TROPHOSPECIES", "WEB"))

write.table(g_df, file="../results/data/g_df.txt", sep="\t", dec=".", quote = F, col.names = T, row.names = F)