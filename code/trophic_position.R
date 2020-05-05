# Determine the trophic position of the species based on predators/prey, 
# including primary producers, herbivores, predators and top predators. 

pwl <- read_delim("../data/kortsch/PairwiseList_seb.txt", delim = " ", col_names = T)


predators <- pwl[,1] %>% 
  distinct()

prey <- pwl[,2] %>%
  distinct()

primary <- pwl %>%
  filter(!PREY %in% PREDATOR) %>%
  distinct(PREY)

both <- pwl %>%
  filter(PREDATOR %in% PREY) %>%
  distinct(PREDATOR)

top <- pwl %>%
  filter(!PREDATOR %in% PREY) %>%
  distinct(PREDATOR)

g_df$TRO_POS <- ""


herb <- g_df %>% filter(TL > 1.1) %>% filter(TL <=2) %>% distinct(ABBREVIATION)

pred <- both %>% filter(!PREDATOR %in% herb$ABBREVIATION) %>%
  distinct(PREDATOR)

for(i in 1:nrow(g_df)){
  if(g_df[i,1] %in% primary$PREY) g_df$TRO_POS[i] <- "primary" 
}


for(i in 1:nrow(g_df)){
  if(g_df[i,1] %in% herb$ABBREVIATION) g_df$TRO_POS[i] <- "herbivore" 
}


for(i in 1:nrow(g_df)){
  if(g_df[i,1] %in% pred$PREDATOR) g_df$TRO_POS[i] <- "predator" 
}


for(i in 1:nrow(g_df)){
  if(g_df[i,1] %in% top$PREDATOR) g_df$TRO_POS[i] <- "top" 
}
write.table(g_df, file="../results/data/g_df.txt", sep="\t", dec=".", quote = F, col.names = T, row.names = F)