#########################################################################################
###                                                                                   ###
###     Script to summarize clusters and open space from stem map and UAS datasets.   ###
###     The script then spatially merges the stem map and UAS datasets at both the    ###
###     tree level and cluster level for subsequent analysis.                         ###
###                                                                                   ###
###     The script is structured to loop through the 11 stem maps and UAS             ###
###     acquisitions in the Black Hills of South Dakota. Outputs from the script      ###
###     include: 1) four panel figures containing diameter and height distributions   ###
###     overlaid for both datasets and stem maps of tree locations for both datasets; ###
###     2) CVS file containing summarized cluster statistics of # of trees, % of      ###
###     stand basal area, height coefficient of variation, and the mean x and y       ###
###     coordinates of the cluster; 3) figure of the distribution of distances from   ###
###     nearest tree overlaid for the two datasets; 4) CSV files of the spatially     ###
###     matched trees evaluated for omission, commission, and true positives, along   ###
###     with precision metrics for UAS height and DBH values; 5) finally CSV files    ###
###     containing spatially matched clusters evaluated for omission, commission, and ###
###     true positives, along with precision metrics of UAS cluster size and height   ###
###     coefficients of variation.                                                    ###
###                                                                                   ###
#########################################################################################

#### Load required packages
library(sf)
library(fpc)       
library(tidyverse)
library(raster)
library(ggplot2)
library(scales)
library(ggpmisc)
library(terra)
library(spatstat.geom)

### Acquire patchwoRk package from Github for patch delineation
### after running this once, patchwoRk can be installed from the library
### devtools::install_github("bi0m3trics/patchwoRk")

library(patchwoRk)

### Set base working directory for analysis
start_time = Sys.time()
rootDir <- "C:/Users/WadeTinkham/Desktop/Black_Hills/point_cloud_processing_BHEF_202306_combined/Cluster_Analysis"
### Create a set of processing folders for outputs ##########################
setwd(rootDir)
dir.create(file.path(rootDir, "/Output"))
dir.create(file.path(rootDir, "/Figures"))

### Minimum Height of overstory trees
min_height = 9

### Load combined tree data, standardize variables, and subset uas and stem map trees
trees = st_read("C:/Users/WadeTinkham/Desktop/Black_Hills/point_cloud_processing_BHEF_202306_combined/final_detected_tree_tops.gpkg")
stands = st_read("C:/Users/WadeTinkham/Desktop/Black_Hills/GIS/Voodoo_Stands.gpkg")

tall.trees <- subset(trees, tree_height_m > min_height)
colnames(stands)[colnames(stands) == "Unit_Numbe"] <- "stand"
  
tall.trees <- st_intersection(tall.trees, stands)
tall.trees$stand <- as.factor(tall.trees$stand)

tree.coordinates <- as.data.frame(as.data.frame(tall.trees$geom)$geometry) %>%
  separate_wider_position(cols = geometry, width = c(2, X = 10, 2, Y = 11), too_many = "drop")
tree.coordinates$X <- as.numeric(tree.coordinates$X)
tree.coordinates$Y <- as.numeric(tree.coordinates$Y)

tall.trees <- cbind(tall.trees, tree.coordinates)

#############################################################################
##### Identify clusters in each stem map plot                           #####
#############################################################################
### Place trees into clusters using and inter-tree distance of 6 m
DBSCAN <- dbscan(tree.coordinates, eps = 6, MinPts = 2)

### append cluster ID to trees
tall.trees$cluster = DBSCAN$cluster
### filter out trees belonging to trees as new dataset
trees.clusters <- subset(tall.trees, cluster > 0)
### Identify trees not belonging to cluster (i.e. isolated individuals)
trees.individuals <- subset(tall.trees, cluster == 0)
### Label individual trees with unique name for that structure
trees.individuals$cluster <- seq(from = max(tall.trees$cluster) + 1, 
                                 to = (max(tall.trees$cluster) + nrow(trees.individuals)),
                                 by = 1)
### combine trees in clusters with individual trees
tall.trees <- rbind(trees.clusters, trees.individuals)

###Summarize clusters
cluster.summary <- data.frame(stand = character(),
                              cluster = numeric(),
                              trees = numeric(),
                              BA_per = numeric(),
                              HT_CV = numeric(),
                              Crown.Area = numeric(),
                              Easting = numeric(),
                              Northing = numeric())

for(i in unique(tall.trees$stand)) {
  stand.trees <- subset(tall.trees, stand == i)
  
  for(j in unique(stand.trees$cluster)) {
    temp.cluster <- subset(stand.trees, cluster == j)
    
    cluster.summary = rbind(cluster.summary,
                             data.frame(stand = i,
                                        cluster = temp.cluster$cluster[1],
                                        trees = nrow(temp.cluster),
                                        BA_per = sum(temp.cluster$basal_area_m2) / sum(stand.trees$basal_area_m2),
                                        HT_CV = sd(temp.cluster$tree_height_m) / mean(temp.cluster$tree_height_m),
                                        Crown.Area = sum(temp.cluster$crown_area_m2),
                                        Easting = mean(temp.cluster$X),
                                        Northing = mean(temp.cluster$Y),
                                        acres = temp.cluster$Ac[1]))
  }
}

### Identify size of each clump
cluster.summary$clump.size <- cut(cluster.summary$trees, 
                                  breaks = c(0,1,4,9,15,10000), 
                                  labels = c("Individual","2-4 trees","5-9 trees","10-15 trees",">15 trees"))

cluster.summary$stand <- as.factor(cluster.summary$stand)
cluster.summary$stand <- factor(cluster.summary$stand, levels=c("1", "2", "3", "4", "9", "10"))
tall.trees$stand <- factor(tall.trees$stand, levels=c("1", "2", "3", "4", "9", "10"))

tall.trees <- tall.trees %>%
  left_join(cluster.summary, by = c("stand", "cluster"))

### Determine nearest neighbor distance for the stand and within clusters
### combination of these metrics has been used to describe relative aggregation
tall.trees$NN_m <- 0
temp.trees <- tall.trees[0,]

for(i in unique(tall.trees$stand)) {
  stand.trees <- subset(tall.trees, stand == i)
  stand.trees$NN_m <- nndist(stand.trees$X, stand.trees$Y, k=1)
  
  temp.trees <- rbind(temp.trees, stand.trees)
}

tall.trees <- temp.trees

### Create boxplots of distance to nearest tree within different clump sizes
ggplot(data = tall.trees, aes(y=NN_m, x = clump.size)) +
  geom_boxplot() +
  labs(x = "Clump Size", 
       y = "Distance to Nearest Neighbor in Clump (m)",
       title = paste0("Histogram of distance to nearest clump neighbor ", "\n" ,min_height," m height threshold")) +
  facet_wrap(vars(stand)) +
  theme_bw()

### Save plots out using stem map plot name to the Figures folder
ggsave(paste0("Figures/Distance_to_nearest_clump_neighbor_", min_height,"m_height_threshold.tiff"), width = 4, height = 3, dpi = 300)

### Save trees that have been assigned to clusters and have metrics of cluster density/variability
write.csv(tall.trees, file = paste0("Output/trees_with_cluster_metrics_", min_height, 
                                    "m_height_threshold.csv"), row.names = FALSE)

### Save summary of cluster density/variability metrics
write.csv(cluster.summary, file = paste0("Output/cluster_metrics_", min_height,
                                         "m_height_threshold.csv"), row.names = FALSE)

### Create map of cluster with unique color for each cluster
number.of.colors = length(unique(tall.trees$clump.size))
cols = rainbow(number.of.colors, s=0.6, v=0.9)[sample(1:number.of.colors,number.of.colors)]

ggplot() +
  geom_sf(data = stands,
          aes(geometry = geom), col = "yellow", fill = "transparent", linewidth = 0.75) +
  geom_sf(data = tall.trees, 
          aes(geometry = geom, col = factor(clump.size)), size = 0.05) +
  coord_sf() +
  scale_color_manual(values = cols) +
  ylim(min(tall.trees$Y - 15), max(tall.trees$Y + 15)) +
  xlim(min(tall.trees$X - 15), max(tall.trees$X + 15)) +
  labs(title = paste0("Clusters using ", min_height," m height threshold"),
       col = "Clump Size") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(paste0("Figures/Cluster Map from ", min_height, " m height threshold.tiff"), width = 8, height = 4, dpi = 300)

### Create stacked bar plot of cluster sizes
ggplot(data = cluster.summary, aes(x = stand)) +
  stat_count(aes(fill = clump.size), position = position_fill()) +
  scale_fill_manual(values = cols) +
  labs(x = "Stand",
       y = "Proportion of Structures",
       fill = "Clump Size",
       title = paste0("Distribution of clump sizes using ", min_height," m height threshold")) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(paste0("Figures/Distribution_of_Clump_Sizes_", min_height,"m_height_threshold.tiff"), width = 6, height = 4)

### Create faceted plot of height histogram for each stand
ggplot(data = tall.trees, aes(x = tree_height_m, group = stand)) +
  facet_wrap(vars(stand), nrow = 2, scales = "free_y") +
  geom_histogram(binwidth = 1, center = 2.5, aes(y=..count../sum(..count..))) +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Height (m)",
       y = "Percent of Trees",
       title = paste0("Height distribution using ", min_height," m height threshold")) +
  xlim(0,30)

ggsave(paste0("Figures/","Height_Histogram_", min_height,"m_height_threshold.tiff"), width = 6, height = 4)


#############################################################################
##### Distance from nearest tree assessment                             #####
#############################################################################

columns <- (round(max(tall.trees$X) + 5,0) - round(min(tall.trees$X) - 5,0)) * 2
rows <- (round(max(tall.trees$Y) + 5,0) - round(min(tall.trees$Y) - 5,0)) * 2

### create a 0.5 m resolution empty raster
r <- raster(ncol = columns,
            nrow = rows,
            crs = st_crs(stands))

### Define raster extent to match buffered UAS trees
extentR_2 <- extent(round(min(tall.trees$X) - 5,0), 
                    round(min(tall.trees$X) - 5,0) + columns, 
                    round(min(tall.trees$Y) - 5,0), 
                    round(min(tall.trees$Y) - 5,0) + rows)
extent(r) <- extent(extentR_2)

### Calculate raster of distance from nearest tree
distance.raster <- distanceFromPoints(r, tree.coordinates)
distance.raster <- mask(distance.raster, mask = stands, updatevalue = NA)

writeRaster(distance.raster, paste0("Output/Distance to tree ", min_height,"m height threshold.tif"), format="GTiff", overwrite=TRUE)

### Create map of distance to nearest tree
distance.raster_spdf <- as(distance.raster, "SpatialPixelsDataFrame")
distance.raster_df <- as.data.frame(distance.raster_spdf)
colnames(distance.raster_df) <- c("value", "x", "y")

ggplot() +
  geom_tile(data = distance.raster_df, aes(x=x, y=y, fill=value)) +
  scale_fill_gradient2(low = "white", mid = "lightyellow", high = "darkgreen", midpoint = 6) +
  geom_sf(data = stands,
          aes(geometry = geom), col = "yellow", fill = "transparent", linewidth = 1) +
  coord_sf() +
  ylim(min(tall.trees$Y) - 10,
       max(tall.trees$Y) + 10) +
  xlim(min(tall.trees$X) - 10,
       max(tall.trees$X) + 10) +
  labs(fill = paste0("Distance From Nearest Tree (m)")) +
  theme_bw() + 
  theme(legend.position = "top")

ggsave(paste0("Figures/","Map of distance to nearest tree", min_height,"m_height_threshold.tiff"), width = 8, height = 4)

### Extract distance from nearest tree based on stand polygons
distance.list <- terra::extract(as(distance.raster, "SpatRaster"), stands)

stands <- c(1,2,3,4,5,6,9,10)
for(i in 1:length(distance.list)) {
  if(length(distance.list[[i]]) > 0) {
    data.frame(dist = distance.list[[i]]) %>%
      ggplot(aes(x=dist)) +
      geom_histogram(binwidth=3, center = 1.5, fill = "darkgrey", col = "black",
                     aes(y=..count../sum(..count..))) +
      scale_y_continuous(labels = percent_format()) +
      labs(x = "Distance from Nearest Tree (m)", 
           y = "Proportion of Stand Area",
           title = paste0("Histogram of distance to nearest tree for ", "\n" ,min_height," m height threshold")) +
      xlim(0, distance.raster@data@max) + 
      theme_bw()
    
    ### Save plots out using stem map plot name to the Figures folder
    ggsave(paste0("Figures/Stand_",stands[i],"_Distance_to_nearest_tree_", min_height,"m_height_threshold.tiff"), width = 4, height = 3, dpi = 300)
  }
}


#############################################################################
##### Gap size assessment                                               #####
#############################################################################
### Reclassify distance raster as gap and not gap
### Based on Matonis and Binkley (2018) who define gaps as areas >6 m from a tree
m <- c(0, 4, NA,  4, Inf, 1)
rclmat <- matrix(m, ncol = 3, byrow = TRUE)
gaps.raster <- reclassify(distance.raster, rclmat)
gaps.raster <- as(gaps.raster, "SpatRaster")
crs(gaps.raster) <- crs(st_read("C:/Users/WadeTinkham/Desktop/Black_Hills/GIS/Voodoo_Stands.gpkg"))

m <- c(0, 9, NA,  9, Inf, 1)
rclmat <- matrix(m, ncol = 3, byrow = TRUE)
openings.raster <- reclassify(distance.raster, rclmat)
openings.raster <- as(openings.raster, "SpatRaster")
crs(openings.raster) <- crs(st_read("C:/Users/WadeTinkham/Desktop/Black_Hills/GIS/Voodoo_Stands.gpkg"))

### Map of Gaps (> 5 m from tree) and openings (> 9 m from tree)
mapview::mapview(gaps.raster, 
                 col.regions= "green",
                 alpha=0.5, na.color = "transparent", 
                 legend = FALSE, 
                 maxpixels=100000000) +
  mapview::mapview(
    tall.trees
    , zcol = "clump.size"
    , cex = 3
    , legend = T
    , label = F
  ) +
  mapview::mapview(openings.raster, 
                   col.regions= "blue",
                   alpha=0.5, na.color = "transparent", 
                   legend = FALSE, 
                   maxpixels=100000000) 

### convert opening raster to polygons
openings.poly <- openings.raster %>%
  patches(directions = 4, allowGaps = FALSE) %>%
  as.polygons(aggregate = TRUE, values = TRUE, na.rm = TRUE) 

crs(openings.poly) <- crs(st_read("C:/Users/WadeTinkham/Desktop/Black_Hills/GIS/Voodoo_Stands.gpkg"))

### Buffer polygons 5 m based on Clyatt et al., 2016
### "... each polygon was buffered by 5 m to expand the polygon to within
### 4 m of neighboring tree piths. The distance of 4 m was used because it 
### approximates the average distance crowns of mature trees in these plots 
### extend into openings."

openings.poly.sf <- sf::st_as_sf(openings.poly)
openings.poly.sf <- sf::st_buffer(openings.poly.sf, dist = 5, nQuadSegs = 30)

#openings.poly <- terra::buffer(openings.poly, width = 5, quadsegs = 90)
# Using terra bufer kept cracshing R

openings.poly <- vect(openings.poly.sf)

### Crop buffered opening polygons by stand boundary
openings.poly <- crop(openings.poly, as(st_read("C:/Users/WadeTinkham/Desktop/Black_Hills/GIS/Voodoo_Stands.gpkg"), "SpatVector"))

### Calculate area attributes
openings.poly$Area_m2 <- expanse(openings.poly, unit = "m")
openings.poly$Area_ha <- openings.poly$Area_m2 / 10000
openings.poly$Area_acre <- openings.poly$Area_ha * 2.47

### Append Stand to gaps polygons
stands = st_read("C:/Users/WadeTinkham/Desktop/Black_Hills/GIS/Voodoo_Stands.gpkg")
colnames(stands)[colnames(stands) == "Unit_Numbe"] <- "stand"

openings.poly <- cbind(openings.poly, extract(vect(st_as_sf(stands)), centroids(openings.poly, inside = TRUE))[,-1])
crs(openings.poly) <- crs(st_read("C:/Users/WadeTinkham/Desktop/Black_Hills/GIS/Voodoo_Stands.gpkg"))

### Assign gap size for each polygon
openings.poly$gap.size <- cut(openings.poly$Area_acre, 
                              breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 10000), 
                              labels = c("> 0.1 ", "0.1-0.2 ", "0.2-0.3 ", "0.3-0.4 ", 
                                         "0.4-0.5 ", "0.5-0.6 ", "0.6-0.7 ", "0.7-0.8 ", 
                                         "0.8-0.9 ", "0.9-1.0 ", "1-2 ", "2-3 ", "> 3 "))

writeVector(openings.poly, paste0("Output/Forest_Openings_", min_height, "m_height_threshold.gpkg"), overwrite = TRUE)

### Create map of forest gaps
ggplot() +
  geom_sf(data = st_as_sf(openings.poly), aes(fill = as.factor(gap.size)), col = "transparent") +
  scale_fill_viridis_d(option = "D") +
  geom_sf(data = stands,
          aes(geometry = geom), col = "yellow", fill = "transparent", linewidth = 0.75) +
  coord_sf() +
  ylim(min(tall.trees$Y) - 10,
       max(tall.trees$Y) + 10) +
  xlim(min(tall.trees$X) - 10,
       max(tall.trees$X) + 10) +
  labs(title = paste0("Gaps when using ", min_height, " m height threshold"), 
       fill = "Gap Size (acres)") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(paste0("Figures/","Map of forest openings ", min_height,"m_height_threshold.tiff"), width = 8, height = 4)

### Distribution of opening sizes
as.data.frame(openings.poly) %>%
  group_by(stand, gap.size) %>%
  summarize(count = n(),
            Gap_ac = sum(Area_acre),
            Stand_ac = mean(Ac)) %>%
  ggplot() +
  geom_col(aes(x = gap.size, y = Gap_ac / Stand_ac, group = stand)) +
  labs(y = "Proportion of Stand", x = "Gap Sizes (acres)") +
  scale_y_continuous(labels = percent_format()) +
  facet_wrap(vars(stand), scales = "free") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 40, vjust = 1, hjust = 1))

ggsave(paste0("Figures/","Relative Gap Size Distribution ", min_height,"m_height_threshold.tiff"), width = 8, height = 4)

as.data.frame(openings.poly) %>%
  group_by(stand, gap.size) %>%
  summarize(count = n(),
            Gap_ac = sum(Area_acre),
            Stand_ac = mean(Ac)) %>%
  ggplot() +
  geom_col(aes(x = gap.size, y = Gap_ac, group = stand)) +
  labs(y = "Acres", x = "Gap Sizes (acres)") +
  facet_wrap(vars(stand), scales = "free") +
  theme_bw() +
  theme(axis.text.x=element_text(angle = 40, vjust = 1, hjust = 1))

ggsave(paste0("Figures/","Cumulative Opening Size Distribution ", min_height,"m_height_threshold.tiff"), width = 8, height = 4)

end_time = Sys.time()

print(paste0("Processing took ", round(end_time - start_time, 2) , " min"))
