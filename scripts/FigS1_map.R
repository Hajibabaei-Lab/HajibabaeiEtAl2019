# Teresita M. Porter, May 31, 2019

library(ggmap) # map
library(ggsn) # map scale bar
library(scales) # comma

# read in sites and coord
siteTable <- read.csv(file="Sites.csv", header=TRUE, sep=",")

siteTable$label <- paste (siteTable$Site, siteTable$Region, sep=" ")

# get map, the higher the zoom, the more detail on the map
map <- get_stamenmap(bbox = c(left = -80.65, bottom = 43.45, 
                              right = -80.50, top = 43.5), 
                              zoom = 13)

p <- ggmap(map) +
  geom_point(data=siteTable, aes(x=Lon, y=Lat), size=1) +
  geom_text(data=siteTable[-6,], aes(x=Lon, y=Lat, label=label), size=4, hjust=-0.1, vjust=-0.1, fontface="bold") +
  geom_text(data=siteTable[6,], aes(x=Lon, y=Lat, label=label), size=4, hjust=1, vjust=1, fontface="bold") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none") +
  scalebar(x.min = -80.65, x.max = -80.5, y.min = 43.45, y.max = 43.5, 
           dist = 1, dist_unit = "km", transform = TRUE, model="WGS84",
           st.bottom = FALSE, st.color = "black", st.size = 3, st.dist = 0.025)

ggsave("FigS1_map.pdf", p, width = 12, height = 4, units = "in")