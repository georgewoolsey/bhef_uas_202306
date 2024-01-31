remove(list = ls()[])
gc()

# setup
setwd("c:/Data/usfs/point_cloud_tree_detection_ex/data/")
delivery_dir = paste0(getwd(),"/point_cloud_processing_BHEF_202306_combined")
# list.files(delivery_dir)

# library
library(tidyverse) # the tidyverse
library(viridis) # viridis colors
library(scales) # work with number and plot scales
# spatial analysis
library(terra) # raster
library(sf) # simple features
# visualization
library(kableExtra)
library(patchwork) # ! better align plots in grid
library(mapview) # interactive html maps
# remove.packages("ggmap")
# devtools::install_github("stadiamaps/ggmap")
library(ggmap)

# option to put satellite imagery as base layer
  mapview::mapviewOptions(
    homebutton = FALSE
    , basemaps = c("Esri.WorldImagery","OpenStreetMap")
  )

# read data
# dtm_rast = terra::rast(paste0(delivery_dir, "/dtm_1m.tif"))
chm_rast = terra::rast(paste0(delivery_dir, "/chm_0.25m.tif"))
dbh_locations_sf = sf::st_read(paste0(delivery_dir, "/bottom_up_detected_stem_locations.gpkg"))
crowns = terra::rast(paste0(delivery_dir, "/top_down_detected_tree_crowns.tif"))
crowns_sf_with_dbh = sf::st_read(paste0(delivery_dir, "/final_detected_crowns.gpkg"))
treetops_sf_with_dbh = sf::st_read(paste0(delivery_dir, "/final_detected_tree_tops.gpkg"))
silv_metrics_temp = readr::read_csv(paste0(delivery_dir, "/final_plot_silv_metrics.csv"))
las_ctg_dta = sf::st_read(paste0(delivery_dir, "/raw_las_ctg_info.gpkg"))

# forest stands shapes
forest_stands = sf::st_read("c:/Data/usfs/bhef_rxfire_plan/data/bhef_harvests.gpkg") %>% 
  dplyr::filter(
    year_id >= year(Sys.time()) - 15 
    & !(treatment_type_grp %in% c("Improvement/Liberation Cut", "Other", "Sanitation Cut"))
  )
# mapview::mapview(forest_stands)
# forest_stands %>% glimpse()
# forest_stands %>% sf::st_drop_geometry() %>% dplyr::count(treatment_type_grp, treatment_type)
##################################################################################
##################################################################################
# summary table
##################################################################################
##################################################################################
  
## summary table
  silv_metrics_temp %>%
    sf::st_drop_geometry() %>%
    dplyr::filter(plotID==silv_metrics_temp$plotID[1]) %>%
    tidyr::pivot_longer(
      cols = -c(plotID)
      , names_to = "measure"
      , values_to = "value"
    ) %>%
    kableExtra::kbl(
      caption = "Silvicultural metrics in both metric and imperial units"
      , digits = 2
    ) %>%
    kableExtra::kable_styling()

#################################################################################
#################################################################################
# Join tree tops with forest stands
#################################################################################
#################################################################################
forest_stands_trees = forest_stands %>%
  dplyr::mutate(
    stand_area_m2 = sf::st_area(.) %>% as.numeric()
    , stand_area_ha = stand_area_m2/10000
  ) %>% 
  sf::st_intersection(las_ctg_dta) %>% 
  dplyr::mutate(
    intrsct_stand_area_m2 = sf::st_area(.) %>% as.numeric()
  ) %>% 
  dplyr::filter(round(intrsct_stand_area_m2, 0) == round(stand_area_m2, 0)) %>% 
  dplyr::select(
    suid, forest_commonname, admin_region_code, activity_name
    , treatment_type, treatment_type_grp, date_compl, year_id, stand_area_m2
  ) %>% 
  sf::st_intersection(treetops_sf_with_dbh)
  
  # forest_stands_trees %>% ggplot() + geom_sf(aes(color=as.factor(suid)), shape=".") + theme(legend.position = "none")
  # forest_stands_trees %>% 
  #   sf::st_drop_geometry() %>% 
  #   dplyr::count(treeID) %>% 
  #   dplyr::filter(n>1) %>% 
  #   summary()

#################################################################################
#################################################################################
# Calculate Silviculture Metrics
#################################################################################
#################################################################################
    # Common silvicultural metrics are calculated for the entire extent.
      # Note, that stand-level summaries can be computed if stand vector data is provided.
      # metrics include:
        # "n_trees"
        # "plot_area_ha"
        # "trees_per_ha"
        # "mean_dbh_cm"
        # "qmd_cm"
        # "mean_tree_height_m"
        # "loreys_height_m"
        # "basal_area_m2"
        # "basal_area_m2_per_ha"

    ### stand-level summaries
      silv_metrics_temp = forest_stands_trees %>%
        sf::st_drop_geometry() %>%
        dplyr::ungroup() %>%
        dplyr::group_by(suid,stand_area_ha) %>%
        dplyr::summarise(
          n_trees = dplyr::n_distinct(treeID)
          , mean_dbh_cm = mean(dbh_cm, na.rm = T)
          , mean_tree_height_m = mean(tree_height_m, na.rm = T)
          , loreys_height_m = sum(basal_area_m2*tree_height_m, na.rm = T) / sum(basal_area_m2, na.rm = T)
          , basal_area_m2 = sum(basal_area_m2, na.rm = T)
          , sum_dbh_cm_sq = sum(dbh_cm^2, na.rm = T)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          trees_per_ha = (n_trees/plot_area_ha)
          , basal_area_m2_per_ha = (basal_area_m2/plot_area_ha)
          , qmd_cm = sqrt(sum_dbh_cm_sq/n_trees)
        ) %>%
        dplyr::select(-c(sum_dbh_cm_sq)) %>%
        # convert to imperial units
        dplyr::mutate(
          dplyr::across(
            .cols = tidyselect::ends_with("_cm")
            , ~ .x * 0.394
            , .names = "{.col}_in"
          )
          , dplyr::across(
            .cols = tidyselect::ends_with("_m")
            , ~ .x * 3.28
            , .names = "{.col}_ft"
          )
          , dplyr::across(
            .cols = tidyselect::ends_with("_m2_per_ha")
            , ~ .x * 4.359
            , .names = "{.col}_ftac"
          )
          , dplyr::across(
            .cols = tidyselect::ends_with("_per_ha") & !tidyselect::ends_with("_m2_per_ha")
            , ~ .x * 0.405
            , .names = "{.col}_ac"
          )
          , dplyr::across(
            .cols = tidyselect::ends_with("_area_ha")
            , ~ .x * 2.471
            , .names = "{.col}_ac"
          )
          , dplyr::across(
            .cols = tidyselect::ends_with("_m2")
            , ~ .x * 10.764
            , .names = "{.col}_ft2"
          )
        ) %>%
        dplyr::rename_with(
          .fn = function(x){dplyr::case_when(
            stringr::str_ends(x,"_cm_in") ~ stringr::str_replace(x,"_cm_in","_in")
            , stringr::str_ends(x,"_m_ft") ~ stringr::str_replace(x,"_m_ft","_ft")
            , stringr::str_ends(x,"_m2_per_ha_ftac") ~ stringr::str_replace(x,"_m2_per_ha_ftac","_ft2_per_ac")
            , stringr::str_ends(x,"_per_ha_ac") ~ stringr::str_replace(x,"_per_ha_ac","_per_ac")
            , stringr::str_ends(x,"_area_ha_ac") ~ stringr::str_replace(x,"_area_ha_ac","_area_ac")
            , stringr::str_ends(x,"_m2_ft2") ~ stringr::str_replace(x,"_m2_ft2","_ft2")
            , TRUE ~ x
          )}
        ) %>%
        dplyr::select(
          "suid"
          , "n_trees"
          , "plot_area_ha"
          , "trees_per_ha"
          , "mean_dbh_cm"
          , "qmd_cm"
          , "mean_tree_height_m"
          , "loreys_height_m"
          , "basal_area_m2"
          , "basal_area_m2_per_ha"
          # imperial
          , "plot_area_ac"
          , "trees_per_ac"
          , "mean_dbh_in"
          , "qmd_in"
          , "mean_tree_height_ft"
          , "loreys_height_ft"
          , "basal_area_ft2"
          , "basal_area_ft2_per_ac"
        )

    ### export tabular
      write.csv(
          silv_metrics_temp
          , paste0(delivery_dir, "/final_plot_silv_metrics.csv")
          , row.names = F
        )

    # this would just be a vector file if available
      silv_metrics_temp = las_ctg_dta %>% 
        dplyr::mutate(
          plotID = "1" # can spatially join to plot vectors if available
        ) %>%
        # join with plot data data
        dplyr::inner_join(
          silv_metrics_temp
          , by = dplyr::join_by("plotID")
        )

    ## summary table
    silv_metrics_temp %>%
      sf::st_drop_geometry() %>%
      dplyr::filter(plotID==silv_metrics_temp$plotID[1]) %>%
      tidyr::pivot_longer(
        cols = -c(plotID)
        , names_to = "measure"
        , values_to = "value"
      ) %>%
      kableExtra::kbl(
        caption = "Silvicultural metrics in both metric and imperial units"
        , digits = 2
      ) %>%
      kableExtra::kable_styling()


##################################################################################
##################################################################################
### Relationship between height and DBH
##################################################################################
##################################################################################
  
### plot
  ggplot(
    data = crowns_sf_with_dbh
    , mapping = aes(y=tree_height_m, x = dbh_cm)
  ) +
  geom_point(
    mapping = aes(color = is_training_data)  
    , alpha = 0.6
    , size = 0.5
  ) + 
  geom_smooth(
    method = "loess"
    , se = F
    , span = 1
    , color = "gray44"
    , alpha = 0.7
  ) +
  labs(
    x = "DBH (cm)"
    , y = "Tree Ht. (m)"
    , color = "Training Data"
    , title = "SfM derived tree height and DBH relationship"
  ) +
  scale_color_manual(values = c("gray", "firebrick")) +
  theme_light() +
  theme(
    legend.position = "bottom"
    , legend.direction = "horizontal"
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = 15, size = 6, alpha = 1))
  )

##################################################################################
##################################################################################
### DBH histogram training vs non-training
##################################################################################
##################################################################################
### plot
  ggplot(
    data = crowns_sf_with_dbh
    , mapping = aes(x = dbh_cm, group = is_training_data, fill = is_training_data)
  ) +
  geom_density(alpha = 0.6, binwidth = 2, color = NA) + 
  labs(
    x = "DBH (cm)"
    , y = "density"
    , fill = "Training Data"
    , title = "SfM derived tree DBH distribution"
  ) +
  scale_fill_manual(values = c("gray", "firebrick")) +
  scale_x_continuous(breaks = scales::extended_breaks(n=20)) +
  theme_light() +
  theme(
    legend.position = "bottom"
    , legend.direction = "horizontal"
  )

##################################################################################
##################################################################################
# chm on html mapview
##################################################################################
##################################################################################
  # aggregate raster and map
  chm_rast %>%
    terra::aggregate(fact=4) %>% 
    stars::st_as_stars() %>% 
    mapview::mapview(
      layer.name = "canopy ht.(m)"
      , alpha.regions = 0.7
    )
##################################################################################
##################################################################################
# chm + trees
##################################################################################
##################################################################################

  ggplot() +
    geom_tile(
      data = chm_rast %>% terra::aggregate(fact=2) %>% as.data.frame(xy=T) %>% dplyr::rename(f=3)
      , mapping = aes(x=x,y=y,fill=f)
    ) +
    geom_sf(
      data = treetops_sf_with_dbh
      # , size = 0.3
      , shape = "."
    ) +
    coord_sf(expand = F) +
    scale_fill_viridis_c(option = "plasma") +
    labs(
      fill="canopy\nht (m)"
      , title = "BHEF Units 1-4 UAS Flights 2023-06"
    ) +  
    theme_light() +
    theme(
      # legend.position =  "top"
      # , legend.direction = "horizontal"
      legend.title = element_text(size = 8, face = "bold")
      , legend.margin = margin(c(0,0,0,0))
      , plot.title = element_text(size = 10, face = "bold") #, hjust = 0.5
      , axis.title = element_blank()
      , axis.text = element_blank()
      , axis.ticks = element_blank()
      , panel.grid = element_blank()
      , plot.margin = margin(0, 0, 0, 0, "cm")
      , plot.caption = element_text(color = "black", hjust = 0, vjust = 3)
    )
  ggplot2::ggsave(
    filename = paste0(delivery_dir,"/plt_chm_treetops.jpeg")
    , plot = ggplot2::last_plot()
    , width = 15
    , height = 13.5
    , units = "in"
    , dpi = "print"
  )
# clean up
  gc()
  remove(list = ls()[grep("_temp",ls())])
##################################################################################
##################################################################################
# try some ggmap !!!!! not working rn
##################################################################################
##################################################################################
if(F){
  ##################hack to align plots for ggmap
  plt_crs = 3857
  ggmap_bbox_fn <- function(map, my_crs=3857) {
      if (!inherits(map, "ggmap")) stop("map must be a ggmap object")
      # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
      # and set the names to what sf::st_bbox expects:
      map_bbox <- setNames(unlist(attr(map, "bb")), c("ymin", "xmin", "ymax", "xmax"))
      # Convert the bbox to an sf polygon, transform it to 3857, 
      # and convert back to a bbox (convoluted, but it works)
      bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), my_crs))
      # Overwrite the bbox of the ggmap object with the transformed coordinates 
      attr(map, "bb")$ll.lat <- bbox_3857["ymin"]
      attr(map, "bb")$ll.lon <- bbox_3857["xmin"]
      attr(map, "bb")$ur.lat <- bbox_3857["ymax"]
      attr(map, "bb")$ur.lon <- bbox_3857["xmax"]
      map
  }
  
  # bounding box
    bb_temp <- 
      # use extent of 
      treetops_sf_with_dbh %>% 
      sf::st_bbox() %>% 
      sf::st_as_sfc() %>% 
      sf::st_transform(crs=5070) %>% 
      sf::st_buffer(as.numeric(100)) %>% 
      sf::st_transform(crs=4326) %>% # same as get_map return
      sf::st_bbox()
      
    # set bbox for get call
    bbox_temp <- c(
      bottom = bb_temp[[2]]
      , top = bb_temp[[4]]
      , right = bb_temp[[3]]
      , left = bb_temp[[1]]
    )
    # get map
    hey_ggmap <- ggmap::get_stadiamap(
      bbox = bbox_temp
      , zoom = 16
      , maptype = "stamen_toner_lite" #"toner-hybrid" # "stamen_terrain" 
      , crop = T
    )
    # ggmap(hey_ggmap)
    # apply align function
    hey_ggmap_aligned <- ggmap_bbox_fn(hey_ggmap, plt_crs) # Use the function
    
    # plot
    # plt_trees = 
      ggmap(hey_ggmap_aligned) +
        geom_sf(
        data= treetops_sf_with_dbh %>% 
          dplyr::filter(tree_height_m>=quantile(crowns_sf_with_dbh$tree_height_m, probs=0.75)) %>% 
          sf::st_transform(crs=4326)
        , size = 0.5
      )
      geom_sf(
        data = crowns_sf_with_dbh %>% 
          dplyr::filter(tree_height_m>=quantile(crowns_sf_with_dbh$tree_height_m, probs=0.75)) %>% 
          sf::st_transform(crs=plt_crs)
        , mapping = aes(fill = tree_height_m)
        , color = NA
      ) 
      +
      geom_sf(
        data= treetops_sf_with_dbh %>% 
          dplyr::filter(tree_height_m>=quantile(crowns_sf_with_dbh$tree_height_m, probs=0.75)) %>% 
          sf::st_transform(crs=plt_crs)
        , size = 0.5
      ) +
      coord_sf(expand = F) +
      scale_fill_viridis_c(option = "plasma") +
      labs(
        fill="canopy\nht (m)"
        , title = "BHEF Units 1-4 UAS Flights 2023-06"
      ) +  
      theme_light() +
      theme(
        # legend.position =  "top"
        # , legend.direction = "horizontal"
        legend.title = element_text(size = 8, face = "bold")
        , legend.margin = margin(c(0,-2,0,0))
        , plot.title = element_text(size = 9, face = "bold") #, hjust = 0.5
        , axis.title = element_blank()
        , axis.text = element_blank()
        , axis.ticks = element_blank()
        , panel.grid = element_blank()
        , plot.margin = margin(0, 0, 0, 0, "cm")
        , plot.caption = element_text(color = "black", hjust = 0, vjust = 3)
      )
    
    plt_trees
    ggplot2::ggsave(
      filename = paste0(delivery_dir,"/plt_chm_treetops.jpeg")
      , plot = plt_trees
      , width = 11
      , height = 10
      , units = "in"
      , dpi = "print"
    )
}
