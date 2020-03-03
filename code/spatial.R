# Test of spatial autocorrelation using the Correlog function from the ncf package
library(rgdal)
library(ncf)

#Loading and converting shape file with coordinates
subshp <- readOGR(dsn="../data/gis/subs.shp") 
subdf <- as_tibble(cbind(as_tibble(subshp@data[,3:6]), subshp@coords)) %>%
  left_join(data_frame(WEB=subnames[-1]), ., by=c("WEB" = "sub"))

#Generate correlogram data based on coordinates and either 
# (normalised) Jaccard distance
# Species overlap
#Setting increments ("bin size") to 50 km
fit_J   <- correlog(x=subdf$coords.x1, y=subdf$coords.x2, z=J_s[-1,-1], latlon = T, increment = 50, resamp = 1000, quiet = T)
fit_ol  <- correlog(x=subdf$coords.x1, y=subdf$coords.x2, z=overlap_matrix[-1,-1], latlon = T, increment = 50, resamp = 1000, quiet = T)

lsize <- 0.2
psize <- 1.5
asize <- 6
ctheme <- theme(axis.title = element_text(size=12), 
                axis.text = element_text(size=10),
                panel.border = element_rect(fill=NA, colour="black"),
                panel.background = element_blank())

c_J <- ggplot(mapping=aes(x=fit_J$mean.of.class, y=as.numeric(fit_J$correlation))) + 
  scale_x_continuous(limits=c(0,1600))+
  geom_point(shape = ifelse(as.numeric(fit_J$p) >= 0.05, 1, 19) , size=psize) +
  geom_hline(yintercept = 0, linetype="dashed") +
  ctheme + 
  labs(x="Distance classes (km)", y="Correlation")
  


c_ol <- ggplot(mapping=aes(x=fit_ol$mean.of.class, y=as.numeric(fit_ol$correlation))) + 
  scale_x_continuous(limits=c(12,1600))+
  geom_point(shape = ifelse(as.numeric(fit_ol$p) >= 0.05, 1, 19) , size=psize) +
  geom_hline(yintercept = 0, linetype="dashed") +
  ctheme + 
  labs(x="Distance classes (km)", y="Correlation")