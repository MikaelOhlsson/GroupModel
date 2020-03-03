library(tidyverse)
#Fix Sebastes spp.
#As it lacks prey species in 14 of the 25 sub regions, it counts as primary producer in those.
#Hence, we add the interactions from all Sebastes species to Sebastes spp.

#######
#Merge the Seb species to SEB_SPP
subr <- read_delim("../data/kortsch/SubregionSpecies.txt", delim = "\t", col_names = F)
subnames <- subr$X1 %>% #List of metaweb + all sub regions
  paste0("sub",.)
subnames[1] <- "Metaweb"
pwl <- read_delim("../data/kortsch/PairwiseList_Kortsch.txt", delim = "\t", col_names = T)
pwl <- pwl[,-1] #Remove column with row numbering 

#Edge list with all the different Sebastes species
sebextra <- pwl %>%
  filter(str_detect(PREDATOR, "SEB_MAR") | str_detect(PREY, "SEB_MAR") | str_detect(PREDATOR, "SEB_MEN") | str_detect(PREY, "SEB_MEN") |  str_detect(PREDATOR, "SEB_VIV") | str_detect(PREY, "SEB_VIV"))

#Renaming
for(i in 1:nrow(sebextra)) {
  if(sebextra$PREDATOR[i] == "SEB_MAR") sebextra[i,1] <- "SEB_SPP"
  if(sebextra$PREY[i] == "SEB_MAR") sebextra[i,2] <- "SEB_SPP"
  if(sebextra$PREDATOR[i] == "SEB_VIV") sebextra[i,1] <- "SEB_SPP"
  if(sebextra$PREY[i] == "SEB_VIV") sebextra[i,2] <- "SEB_SPP"
  if(sebextra$PREDATOR[i] == "SEB_MEN") sebextra[i,1] <- "SEB_SPP"
  if(sebextra$PREY[i] == "SEB_MEN") sebextra[i,2] <- "SEB_SPP"
}

pwx <- rbind(pwl, sebextra)  #Bind SEB_SPP with extra interactions to the full interaction list (which includes previous SEB_SPP ints)
 
pwx <- unique(pwx[,1:2])  #Remove duplicate interactions
  
write.table(pwx, file = "../data/kortsch/PairwiseList_seb.txt", sep = " ", row.names = F, quote = F)
