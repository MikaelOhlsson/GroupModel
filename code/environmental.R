library(vegan) #adonis and vegdist
# Running Permanovas to test group structure vs environmental variables

perm <- 10000 #Permutations to use in permanonvas

# # # Permanovas for testing environmental variables # # #
env_anova <- adonis(J_s[-1,-1]~subdf$wctemp+subdf$depth, permutations = perm, method="euclidean") 
