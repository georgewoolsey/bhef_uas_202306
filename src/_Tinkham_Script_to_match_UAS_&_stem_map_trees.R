#########################################################################################
###                                                                                   ###
###     This script matches stem map trees with segmented trees from UAS remote       ###
###     sensing products. It then calculates both individual tree precision metrics   ###
###     and overall summary metrics by canopy stratum layers.                         ###
###                                                                                   ###
#########################################################################################

library(pacman)
pacman::p_load(sp, tidyverse, rgeos, gridExtra, ggpubr, ggbreak, viridis)

### Define root directory and create subsequent output directories
rootDir <- "C:/Users/WadeTinkham/Documents/Research/UAS_Projects/DBH_Modeling_Paper"
setwd(rootDir)

dir.create(file.path(rootDir, "/Analysis"))
dir.create(file.path(rootDir, "/Analysis/Outputs"))
dir.create(file.path(rootDir, "/Analysis/Figures"))

### Read in the stem mapped tree data
map_trees <- read.csv("Stem_Map_Trees.csv")

### Read in UAS detected trees
uas_trees.original <- read.csv("final_uas_trees_unfiltered.csv")
uas_trees.original <- read.csv("final_uas_trees_filtered.csv")


### Variable standardization (This will need to customized for your data)
uas_trees <- uas_trees.original

uas_trees$ht.uas <- uas_trees$tree_heigh
uas_trees$treeID.UAS <- uas_trees$treeID
uas_trees$x.uas <- uas_trees$X
uas_trees$y.uas <- uas_trees$Y
uas_trees$uas.dbh <- uas_trees$dbh * 100  ### only use for unfiltered uas trees

#UAS_trees$treeID.UAS <- seq(1, nrow(uas_trees.original), 1)
uas_trees <- uas_trees %>% mutate(uas.dbh = replace_na(uas.dbh, 0))
is.na(uas_trees$uas.dbh) <- uas_trees$uas.dbh >= 1000


stem_map_trees <- data.frame(treeID.map = map_trees$Tree, 
                             y.map = map_trees$Northing, 
                             x.map = map_trees$Easting, 
                             Species = map_trees$Species,
                             dbh.map = map_trees$DBH_cm,     
                             cbh.map = map_trees$CBHc_18 / 3.28,          ### I assumed CBH was in feet and needed converted to meters  
                             ht.map = map_trees$HT_m)

###############################################################################
### The below set of parameters and pathways need to be defined prior       ###
### to running.                                                             ###
###############################################################################

### True tree file to only include trees above a specific height
min_tree_height = 1.37

### Find all nearest neighbors within maximum distance and then filter them based on 
### distance and height error 
max_dist = 4
max_height_error = 3

#######################################################################################
###     You may be able to skip this if your data is clean                          ###
#######################################################################################

#uas_trees <- subset(uas_trees, ht.uas > (min_tree_height - max_height_error))
uas_trees <- uas_trees[order(uas_trees$treeID.UAS),]
stem_map_trees <- subset(stem_map_trees, ht.map > min_tree_height)
stem_map_trees <- subset(stem_map_trees, dbh.map > 0)

#######################################################################################
###     Everything past here is automated and can be run once the above             ###
###     parameters are taken care of.                                               ###
#######################################################################################

### Convert stems to spatial object for identifying candidate neighbors
target_map <- SpatialPointsDataFrame(
  coords = cbind(stem_map_trees$x.map, stem_map_trees$y.map),
  data = stem_map_trees)

### Setup output dataframe for tree matches
matches <- data.frame(treeID.map = numeric(), y.map = numeric(), x.map = numeric(), dbh.map = numeric(), ht.map = numeric(), 
                      treeID.UAS = numeric(), x.uas = numeric(), y.uas = numeric(), dist = numeric(), 
                      ht.uas = numeric(), uas.dbh = numeric(), ht_error = numeric(), ht_abs_error = numeric(), 
                      dbh_error = numeric(), dbh_abs_error = numeric(), Summary = character())

###############################################################################
### Loop through comparing the individual extraction trees and stems        ###
### before matching extracted trees and stem map trees, which then have     ###
### summary statistics appended to the output database.                     ###
###############################################################################  
### Convert to spatial object for identifying candidate neighbors
target_trees <- SpatialPointsDataFrame(
  coords = cbind(stem_map_trees$x.map, stem_map_trees$y.map),
  data = stem_map_trees)

for(k in uas_trees$treeID.UAS){
  ### Loop through identifying the best stem map match for each uas tree, once a matching stem 
  ### map tree is found, it is removed from the stem map tree so that it cannot be selected 
  ### again. This process ensures that the best one-to-one match is found for each uas tree.
  tmp_plot <- subset(uas_trees, treeID.UAS == k)
  tmp_plot <- SpatialPointsDataFrame(
    coords = cbind(tmp_plot$x.uas, tmp_plot$y.uas),
    data = tmp_plot)
  tmp_plot_buff <- gBuffer(tmp_plot, width = max_dist)
  tmp_trees <- target_trees[tmp_plot_buff,]
  
  if(dim(tmp_trees)[1] == 0){
    stem_map_trees <- as.data.frame(target_trees)
    stem_map_trees <- stem_map_trees[,c(1:7)]
    target_trees <- SpatialPointsDataFrame(
      coords = cbind(stem_map_trees$x.map, stem_map_trees$y.map),
      data = stem_map_trees)
    
    print(k)
  } else {
    tmp_trees$x.uas = tmp_plot$x.uas
    tmp_trees$y.uas = tmp_plot$y.uas
    tmp_trees$diffX = tmp_trees@data$x.uas - tmp_trees@coords[,1]
    tmp_trees$diffY = tmp_trees@data$y.uas - tmp_trees@coords[,2]
    tmp_trees$dist = with(tmp_trees@data, sqrt(diffX^2 + diffY^2))
    tmp_trees$treeID.UAS = k
    tmp_trees$ht.uas = tmp_plot@data$ht.uas
    tmp_trees$uas.dbh = tmp_plot@data$uas.dbh
    tmp_trees$ht_error = tmp_trees$ht.uas - tmp_trees$ht.map
    tmp_trees$ht_abs_error = abs(tmp_trees$ht_error)
    tmp_trees$dbh_error = tmp_trees$uas.dbh - tmp_trees$dbh.map
    tmp_trees$dbh_abs_error = abs(tmp_trees$dbh_error)

    tmp_trees <- as.data.frame(tmp_trees@data)
    
    ### Sort to find the smallest height error for each set of matches and merge all uas 
    ### and stem map variables to a single temporary file
    tmp_trees <- tmp_trees[order(tmp_trees$ht_abs_error),]
    tmp_plot <- as.data.frame(tmp_plot)
    
    if(tmp_trees$ht_abs_error[1] <= max_height_error){
      ### Keep closest stem map tree to uas tree and bind with other matches
      tmp_dat <- data.frame(tmp_trees[1,], Summary = "True Positive")
      matches <- rbind(matches, tmp_dat)
      
      ### If a good stem map tree is identified as a match, call it out as the tree_target 
      ### and use this to remove that tree from the stem map dataframe and then reconvert 
      ### this to a spatial object.
      tree_target <- tmp_dat$treeID.map[1]
      stem_map_trees <- as.data.frame(target_trees)
      stem_map_trees <- stem_map_trees[,c(1:7)]
      
      stem_map_trees <- subset(stem_map_trees, !(stem_map_trees$treeID.map == tree_target))
      
      target_trees <- SpatialPointsDataFrame(
        coords = cbind(stem_map_trees$x.map, stem_map_trees$y.map),
        data = stem_map_trees)
    } else {
      stem_map_trees <- as.data.frame(target_trees)
      stem_map_trees <- stem_map_trees[,c(1:7)]
      
      target_trees <- SpatialPointsDataFrame(
        coords = cbind(stem_map_trees$x.map, stem_map_trees$y.map),
        data = stem_map_trees)
    }
  }
  print(k)
}


### Append all unmatched UAS trees to new dataframe as commissions
commissions = data.frame(filter(uas_trees, !(treeID.UAS %in% matches$treeID.UAS)),
                         Summary = "Commission",
                         treeID.map = NA, y.map = NA, x.map = NA, dbh.map = NA, ht.map = NA, dist = NA, 
                         diffX = NA, diffY = NA, Species = NA, cbh.map = NA,
                         ht_error = NA, ht_abs_error = NA, dbh_error = NA, dbh_abs_error = NA)
  
matches <- rbind(matches, commissions[,c(62:67,68:81)])  ### for unfiltered uas trees
matches <- rbind(matches, commissions[,c(63,65:69,70:83)])  ### for filtered uas trees

#matches <- rbind(matches, commissions[,c(1:4,56,61:67,72:75)])

### Append all unmatched field stem mapped trees to dataframe as omissions
omissions = data.frame(filter(stem_map_trees, !(treeID.map %in% matches$treeID.map)),
                       Summary = "Omission", dist = NA, x.uas = NA, y.uas = NA,
                        diffX = NA, diffY = NA,
                       treeID.UAS = NA, ht.uas = NA, uas.dbh = NA, ht_error = NA, ht_abs_error = NA, 
                       dbh_error = NA, dbh_abs_error = NA)
#omissions <- omissions[,c(1:3,5,7:18)]

matches <- rbind(matches, omissions)

matches$uas.dbh[matches$uas.dbh == 0] <- NA
matches$dbh_error[is.na(matches$uas.dbh)] <- NA
matches$dbh_abs_error[is.na(matches$uas.dbh)] <- NA


### Calculate summary statistics and append to output dataframe 
### Setup output dataframe for summary statistics
size.class.output <- data.frame(dataset = character(), stratum = character(), map_count = numeric(), 
                                uas_count = numeric(), uas_dbh_n = numeric(), ME_ht = numeric(), 
                                RMSE_ht = numeric(), RMSEP_ht = numeric(), TP = numeric(), 
                                FN = numeric(), FP = numeric(), ME_dbh = numeric(), 
                                RMSE_dbh = numeric(), RMSEP_dbh = numeric(),
                                ME_dbh_predict = numeric(), RMSE_dbh_predict = numeric(),
                                RMSEP_dbh_predict = numeric())
  
output <- data.frame(map_count = nrow(matches) - length(matches$treeID.map[matches$Summary == "Commission"]),
                     uas_count = length(uas_trees$ht.uas[uas_trees$ht.uas > min_tree_height]),
                     uas_dbh_n = nrow(subset(matches, uas.dbh > 0)),
                     ME_ht = mean(matches$ht_error[matches$Summary == "True Positive"], na.rm=TRUE), 
                     RMSE_ht = mean(matches$ht_error[matches$Summary == "True Positive"]^2)^0.5, 
                     RMSEP_ht = mean((matches$ht_error[matches$Summary == "True Positive"] / matches$ht.map[matches$Summary == "True Positive"] * 100)^2, na.rm=TRUE)^0.5, 
                     ME_dbh = mean(matches$dbh_error[matches$Summary == "True Positive"], na.rm=TRUE), 
                     RMSE_dbh = mean(matches$dbh_error[matches$Summary == "True Positive"]^2, na.rm=TRUE)^0.5, 
                     RMSEP_dbh = mean((matches$dbh_error[matches$Summary == "True Positive"] / matches$dbh.map[matches$Summary == "True Positive"] * 100)^2, na.rm=TRUE)^0.5, 
                     TP = length(matches$treeID.map[matches$Summary == "True Positive"]) / (length(matches$treeID.map[matches$Summary == "True Positive"]) + length(matches$treeID.map[matches$Summary == "Omission"])),
                     FN = length(matches$treeID.map[matches$Summary == "Omission"]) / (length(matches$treeID.map[matches$Summary == "True Positive"]) + length(matches$treeID.map[matches$Summary == "Omission"])), 
                     FP = length(matches$treeID.UAS[matches$Summary == "Commission"]) / (length(matches$treeID.UAS[matches$Summary == "Commission"]) + length(matches$treeID.map[matches$Summary == "True Positive"])))
  
for(m in seq(5, round(max(uas_trees$ht.uas), 0)+1, 5)){
  matches.temp <- matches
  matches.temp$ht.temp <- pmax(matches.temp$ht.map, matches.temp$ht.uas, na.rm = TRUE)
  matches.subset <- subset(matches.temp, ht.temp >= m-5 & ht.temp < m)
  UAS_trees.subset <- subset(uas_trees, ht.uas >= m-4.4 & ht.uas < m+0.6)
    
  size.class.output <- rbind(size.class.output,
                             data.frame(class.upper.bound = m,
                                        map_count = nrow(matches.subset) - length(matches.subset$treeID.map[matches.subset$Summary == "Commission"]),
                                        uas_count = length(UAS_trees.subset$ht.uas[UAS_trees.subset$ht.uas > min_tree_height]),
                                        uas_dbh_n = length(matches.subset$uas.dbh) - sum(is.na(matches.subset$uas.dbh)),
                                        ME_ht = mean(matches.subset$ht_error[matches.subset$Summary == "True Positive"], na.rm=TRUE), 
                                        RMSE_ht = mean((matches.subset$ht.map[matches.subset$Summary == "True Positive"] - matches.subset$ht.uas[matches.subset$Summary == "True Positive"])^2, na.rm=TRUE) ^0.5, 
                                        RMSEP_ht = mean(((matches.subset$ht.map[matches.subset$Summary == "True Positive"] - matches.subset$ht.uas[matches.subset$Summary == "True Positive"]) / matches.subset$ht.map[matches.subset$Summary == "True Positive"] * 100)^2, na.rm=TRUE) ^0.5, 
                                        ME_dbh = mean(matches.subset$dbh_error[matches.subset$Summary == "True Positive"], na.rm=TRUE),                                           RMSE_dbh = mean((matches.subset$dbh.map[matches.subset$Summary == "True Positive"] - matches.subset$uas.dbh[matches.subset$Summary == "True Positive"])^2, na.rm=TRUE) ^0.5, 
                                        RMSEP_dbh = mean(((matches.subset$dbh.map[matches.subset$Summary == "True Positive"] - matches.subset$uas.dbh[matches.subset$Summary == "True Positive"]) / matches.subset$dbh.map[matches.subset$Summary == "True Positive"] * 100)^2, na.rm=TRUE) ^0.5, 
                                        TP = length(matches.subset$treeID.map[matches.subset$Summary == "True Positive"]) / (length(matches.subset$treeID.map[matches.subset$Summary == "True Positive"]) + length(matches.subset$treeID.map[matches.subset$Summary == "Omission"])),
                                        FN = length(matches.subset$treeID.map[matches.subset$Summary == "Omission"]) / (length(matches.subset$treeID.map[matches.subset$Summary == "True Positive"]) + length(matches.subset$treeID.map[matches.subset$Summary == "Omission"])), 
                                        FP = length(matches.subset$treeID.UAS[matches.subset$Summary == "Commission"]) / (length(matches.subset$treeID.UAS[matches.subset$Summary == "Commission"]) + length(matches.subset$treeID.map[matches.subset$Summary == "True Positive"]))))
}
  
output <- rbind(data.frame(output, class.upper.bound = "Overall"), 
                size.class.output)
  
height.model <- lm(data = matches, ht.map ~ ht.uas)
height.model.summary <- summary(height.model)
height.model.lab <- paste0("Observed = ", round(height.model$coefficients[1],3), " (",round(height.model.summary$coefficients[1,2],2),") + Predicted x ", round(height.model$coefficients[2],3), " (",round(height.model.summary$coefficients[2,2],3),")")
height.r.labs <- paste0('R^2 = ', round(summary(height.model)$adj.r.squared,4),'    RMSE = ', round(output$RMSE_ht[output$class.upper.bound == "Overall"],2), " m")
  
height.error <- ggplot(matches, aes(x = ht.uas, y = ht.map)) +
  geom_point(alpha = 0.5)+
  labs(x = "Stem Map Heights (m)", y = "UAS Detected Heights (m)") +
  geom_smooth(method = 'lm', formula = y~x, color = "black", linetype = 2) +
  xlim(0, 28) +
  ylim(0, 28) +
  annotate("text", x = 7.2, y = 25.8, label = height.r.labs, size = 2) +
  annotate("text", x = 11.4, y = 27, label = height.model.lab, size = 2)
height.error
  
ggsave("Analysis/Figures/unfiltered_height_comparison.jpg", width = 3, height = 3, dpi = 300)
ggsave("Analysis/Figures/filtered_height_comparison.jpg", width = 3, height = 3, dpi = 300)


dbh.model <- lm(data = subset(matches, matches$uas.dbh > 0), dbh.map ~ uas.dbh)
dbh.model.summary <- summary(dbh.model)
dbh.model.lab <- paste0("Observed = ", round(dbh.model$coefficients[1],3), " (",round(dbh.model.summary$coefficients[1,2],2),") + Predicted x ", round(dbh.model$coefficients[2],3), " (",round(dbh.model.summary$coefficients[2,2],2),")")
dbh.r.labs <- paste0('R^2 = ', round(dbh.model.summary$adj.r.squared,4), '    RMSE = ', round(output$RMSE_dbh[output$class.upper.bound == "Overall"],1)," cm")
  
dbh.error <- ggplot(matches, aes(x = uas.dbh, y = dbh.map)) +
  geom_point(alpha = 0.5)+
  labs(y = "Stem Map DBH (cm)", x = "UAS Detected DBH (cm)") +
  geom_smooth(method = 'lm', formula = y~x, color = "black", linetype = 2) +
  xlim(0, 70) +
  ylim(0, 70) +
  annotate("text", x = 20, y = 66.6, label = dbh.r.labs, size = 2) +
  annotate("text", x = 20.6, y = 69, label = dbh.model.lab, size = 2)
dbh.error  

ggsave(plot=dbh.error, "Analysis/Figures/unfiltered_dbh_comparison.tiff", width = 3, height = 3, dpi = 300)
ggsave(plot=dbh.error, "Analysis/Figures/filtered_dbh_comparison.tiff", width = 3, height = 3, dpi = 300)


output$FP[output$FP<0] <- 0
output[is.na(output)] = 0

### Calculate final statistics based on True Positive, False Negative, and False Positive
output$recall = output$TP / (output$TP + output$FN)
output$precision = output$TP / (output$TP + output$FP) 
output$f_score = 2 * output$recall * output$precision / (output$recall + output$precision)

### Standardize factor levels and scale extraction rates to be percent
output$class.upper.bound <- plyr::revalue(output$class.upper.bound, c("5"="1.37 - 5", "10"="5 - 10", "15"="10 - 15",
                                                                              "20"="15 - 20", "25"="20 - 25", "30"="25 - 30",
                                                                              "Overall"="Overall"))

output$class.upper.bound <- factor(output$class.upper.bound, 
                                       levels = c("1.37 - 5", "5 - 10", "10 - 15",
                                                  "15 - 20", "20 - 25", "> 25",
                                                  "Overall"))

output <- droplevels(output[!output$class.upper.bound == '> 25',])
output$FP <- output$FP * 100
output$TP <- output$TP * 100
output$FN <- output$FN * 100

#########################################################################################
### Write out csv files with the first containing summary precision metrics related   ###
### to tree extraction and tree metrics and the second containing all True Positive,  ###
### False Positive, and False Negative tree matches. The UAS extracted trees are      ###
### represented by TP and FP, while stem map trees are the TP and FN.                 ###
#########################################################################################
write.csv(output, "Analysis/Outputs/unfiltered_uas_summary_precision_metrics.csv", row.names = FALSE)
write.csv(matches, "Analysis/Outputs/unfiltered_combined_uas_&_stem_map_trees.csv", row.names = FALSE)

write.csv(output, "Analysis/Outputs/filtered_uas_summary_precision_metrics.csv", row.names = FALSE)
write.csv(matches, "Analysis/Outputs/filtered_combined_uas_&_stem_map_trees.csv", row.names = FALSE)



#########################################################################################
### Generate summary graphics of tree metric accuracy and precision and then of tree  ###
### extraction confusion matrix with F-score.                                         ###
#########################################################################################

### Plotting Height and DBH Errors
ht.me.plot <- ggplot(data = output, aes(x = class.upper.bound, y = ME_ht)) +
  geom_col() +
  geom_hline(yintercept = 0) +
  ylab("Height Mean Error (m)") +
  xlab("Height Size Classes (m)") +
  ylim(-1,1) +
  theme_bw()
ht.me.plot

ht.rmse.plot <- ggplot(data = output, aes(x = class.upper.bound, y = RMSE_ht)) +
  geom_col() +
  ylab("Height RMSE (m)") +
  xlab("Height Size Classes (m)") +
  #ylim(0,2) +
  theme_bw()
ht.rmse.plot

dbh.me.plot <- ggplot(data = output, aes(x = class.upper.bound, y = ME_dbh)) +
  geom_col() +
  geom_hline(yintercept = 0) +
  ylab("DBH Mean Error (cm)") +
  xlab("Height Size Classes (m)") +
  #ylim(-1, 5) +
  theme_bw()
dbh.me.plot

dbh.rmse.plot <- ggplot(data = output, aes(x = class.upper.bound, y = RMSE_dbh)) +
  geom_col() +
  ylab("DBH RMSE (cm)") +
  xlab("Height Size Classes (m)") +
  #ylim(0, 10) +
  theme_bw()
dbh.rmse.plot

combined.plot <- grid.arrange(ht.me.plot,
                              ht.rmse.plot,
                              dbh.me.plot,
                              dbh.rmse.plot,
                              nrow=4, ncol=1)
as_ggplot(combined.plot)

ggsave("Analysis/Figures/unfiltered_height_&_DBH_errors_by_class.jpg", width = 4, height = 8, dpi = 300)
ggsave("Analysis/Figures/filtered_height_&_DBH_errors_by_class.jpg", width = 4, height = 8, dpi = 300)


### Tree Extraction Precision
tp.score.plot <- ggplot(data = output, aes(x = class.upper.bound, y = TP)) +
  geom_col() +
  ylab("True Positive (%)") +
  xlab("Height Size Classes (m)") +
  ylim(0,100) +
  theme_bw()
tp.score.plot

fp.score.plot <- ggplot(data = output, aes(x = class.upper.bound, y = FP)) +
  geom_col() +
  ylab("False Positive (%)") +
  xlab("Height Size Classes (m)") +
  ylim(0,100) +
  theme_bw()
fp.score.plot

fn.score.plot <- ggplot(data = output, aes(x = class.upper.bound, y = FN)) +
  geom_col() +
  ylab("False Negative (%)") +
  xlab("Height Size Classes (m)") +
  ylim(0,100) +
  theme_bw()
fn.score.plot

f.score.plot <- ggplot(data = output, aes(x = class.upper.bound, y = f_score)) +
  geom_col() +
  ylab("F-Score") +
  xlab("Height Size Classes (m)") +
  ylim(0,1) +
  theme_bw()
f.score.plot

combined.plot <- grid.arrange(f.score.plot,
                              tp.score.plot,
                              fp.score.plot,
                              fn.score.plot,
                              nrow=4, ncol=1)
as_ggplot(combined.plot)

ggsave("Analysis/Figures/Extraction_Quality_by_Height_Class.jpg", width = 4, height = 8, dpi = 300)



### Generate plots of UAS extracted DBH distribution
unfiltered_matches <- read.csv("Analysis/Outputs/unfiltered_combined_uas_&_stem_map_trees.csv")
filtered_matches <- read.csv("Analysis/Outputs/filtered_combined_uas_&_stem_map_trees.csv")

dbh.data <- data.frame(Data_Source = c(rep("Stem Map", nrow(subset(unfiltered_matches, unfiltered_matches$dbh.map > 0))),
                                       rep("UAV Unfiltered", nrow(subset(unfiltered_matches, unfiltered_matches$uas.dbh > 0))),
                                       rep("UAV Filtered", nrow(subset(filtered_matches, filtered_matches$uas.dbh > 0)))),
                       DBH = c(unfiltered_matches$dbh.map[is.finite(unfiltered_matches$dbh.map)],
                               unfiltered_matches$uas.dbh[is.finite(unfiltered_matches$uas.dbh)],
                               filtered_matches$uas.dbh[is.finite(filtered_matches$uas.dbh)]))
dbh.data$Data_Source <- as.factor(dbh.data$Data_Source)
dbh.data$Data_Source <- factor(dbh.data$Data_Source, levels = c("Stem Map", "UAV Unfiltered", "UAV Filtered"))

map_distributions <- ggplot(data = dbh.data, aes(x = DBH, after_stat(count/2.7), group = Data_Source, fill = Data_Source)) +
  geom_histogram(binwidth = 2.5, center = 1.25, alpha = 0.4, position = 'identity', color = "grey16") +
 # geom_histogram(data = matches, aes(x = uas.dbh, after_stat(count/4)), color = "blue", fill = "lightblue", binwidth = 2.5, alpha=0.5, center = 1.25) +
  labs(y = "Trees per Hectare", x = "DBH (cm)", fill = "") +
  #ylim(0,38) +
  #xlim(0,346) +
  xlim(0,160) +
  scale_y_break(c(38,63), scales = 0.5) +
  scale_fill_viridis(discrete = TRUE) +
  #scale_x_break(c(160,230), scales = 0.5) +
  #scale_colour_manual(name = "Data Source", values = c("grey", "black"), limits = c("Stem Map", "UAS Extrated")) +
  #scale_fill_manual(name = "Data Source", values = c("#e5f5e0", "#a1d99b", "#31a354"), limits = c("Stem Map", "UAV Unfiltered", "UAV Filtered")) +
  #guides(colour = guide_legend(override.aes = list(fill = c("lightgrey", "grey50")))) +
  theme_bw() +
  theme(legend.position = "top", legend.title = element_text(size=10), legend.text = element_text(size=8)) 
  
map_distributions
  
ggsave("Analysis/Figures/Extracted_DBH_Distribution_yrange.jpg", width = 6, height = 4, dpi = 300)
ggsave("Analysis/Figures/Extracted_DBH_Distribution_xrange.jpg", width = 6, height = 2.9, dpi = 300)

### Extraction Proportions

extraction_summary <- data.frame(Data_Source = c(rep("Stem Map", 3), rep("UAS Unfiltered", 3), rep("UAS Filtered", 3)),
                                 Size_Class = rep(c("< 20 cm", "20 - 40 cm", "> 40 cm"), 3),
                                 Count = c(nrow(map_trees[map_trees$DBH_cm < 20,]), 
                                           nrow(map_trees[map_trees$DBH_cm > 20 & map_trees$DBH_cm < 40,]), 
                                           nrow(map_trees[map_trees$DBH_cm > 40,]), 
                                           (nrow(unfiltered_matches[unfiltered_matches$uas.dbh < 20 & unfiltered_matches$uas.dbh > 0.1,]) - nrow(unfiltered_matches[unfiltered_matches$uas.dbh == "NA",])),
                                           (nrow(unfiltered_matches[unfiltered_matches$uas.dbh > 20 & unfiltered_matches$uas.dbh < 40,]) - nrow(unfiltered_matches[unfiltered_matches$uas.dbh == "NA",])),
                                           (nrow(unfiltered_matches[unfiltered_matches$uas.dbh > 40,]) - nrow(unfiltered_matches[unfiltered_matches$uas.dbh == "NA",])),
                                           (nrow(filtered_matches[filtered_matches$uas.dbh < 20 & filtered_matches$uas.dbh > 0.1,]) - nrow(filtered_matches[filtered_matches$uas.dbh == "NA",])),
                                           (nrow(filtered_matches[filtered_matches$uas.dbh > 20 & filtered_matches$uas.dbh < 40,]) - nrow(filtered_matches[filtered_matches$uas.dbh == "NA",])),
                                           (nrow(filtered_matches[filtered_matches$uas.dbh > 40,]) - nrow(filtered_matches[filtered_matches$uas.dbh == "NA",]))),
                                 Percentage = c(100,100,100,
                                                (nrow(unfiltered_matches[unfiltered_matches$uas.dbh < 20 & unfiltered_matches$uas.dbh > 0.1,]) - nrow(unfiltered_matches[unfiltered_matches$uas.dbh == "NA",])) / nrow(map_trees[map_trees$DBH_cm < 20,]) * 100,
                                                (nrow(unfiltered_matches[unfiltered_matches$uas.dbh > 20 & unfiltered_matches$uas.dbh < 40,]) - nrow(unfiltered_matches[unfiltered_matches$uas.dbh == "NA",]))/ nrow(map_trees[map_trees$DBH_cm > 20 & map_trees$DBH_cm < 40,]) * 100,
                                                (nrow(unfiltered_matches[unfiltered_matches$uas.dbh > 40,]) - nrow(unfiltered_matches[unfiltered_matches$uas.dbh == "NA",])) / nrow(map_trees[map_trees$DBH_cm > 40,]) * 100,
                                                (nrow(filtered_matches[filtered_matches$uas.dbh < 20 & filtered_matches$uas.dbh > 0.1,]) - nrow(filtered_matches[filtered_matches$uas.dbh == "NA",])) / nrow(map_trees[map_trees$DBH_cm < 20,]) * 100,
                                                (nrow(filtered_matches[filtered_matches$uas.dbh > 20 & filtered_matches$uas.dbh < 40,]) - nrow(filtered_matches[filtered_matches$uas.dbh == "NA",])) / nrow(map_trees[map_trees$DBH_cm > 20 & map_trees$DBH_cm < 40,]) * 100,
                                                (nrow(filtered_matches[filtered_matches$uas.dbh > 40,]) - nrow(filtered_matches[filtered_matches$uas.dbh == "NA",])) / nrow(map_trees[map_trees$DBH_cm > 40,]) * 100
                                                ))

write.csv(extraction_summary, "Analysis/Outputs/UAS_Extraction_Summary.csv", row.names = FALSE)

