#=============================================================================
# scripts/simulate_progression.R
# Grounded 30-Year Spatiotemporal Network Outbreak Projection & Visualization
# Production-Grade Mathematically Aligned Edition
#=============================================================================

library(tidyverse)
library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(gganimate)
library(ggspatial) 
library(av)
library(ggnewscale)

#-----------------------------------------------------------------------------
# 1. Ingest Pipeline Spatial Assets & Calibrated Parameters ----
#-----------------------------------------------------------------------------
cat("-> Verifying structural pipeline file dependencies...\n")
if(!file.exists("outputs/nodes.rds") || !file.exists("outputs/movement_matrix.rds")) {
  stop("CRITICAL PIPELINE FAULT: Network nodes or adjacency matrices are missing. Run previous scripts first.")
}
if(!file.exists("outputs/optimal_abc_parameters.rds")) {
  stop("CRITICAL PIPELINE FAULT: Parameters uncalibrated. Run 'scripts/ibm_smc_abc.R' first.")
}

# Load the empirical network nodes and movement configurations
empirical_nodes           <- readRDS("outputs/nodes.rds")
empirical_movement_matrix <- readRDS("outputs/movement_matrix.rds")
env_matrix                <- unname(as.matrix(empirical_movement_matrix))

# Load ABC optimized parameter metrics matching the SMC estimation outputs
abc_params             <- readRDS("outputs/optimal_abc_parameters.rds")
transmission_scale_val <- abc_params %>% 
  filter(parameter == "transmission_scale") %>% 
  pull(optimal_value)

if(length(transmission_scale_val) == 0) {
  stop("PIPELINE MISMATCH: 'transmission_scale' target parameter missing from calibration file.")
}

N_nodes     <- nrow(empirical_nodes)
total_years <- 30

# ============================================================================
# TARGET FIX: Sourcing & Realignment for Perfect Border Edge Snapping
# ============================================================================
cat(">>> Ingesting pre-cropped landscape boundaries from cassMap.R...\n")
source("scripts/cassMap.R") 

cat("-> Re-projecting vector assets to UTM 30N (EPSG:32630)...\n")
country_borders_utm <- st_transform(three_countries_proj, 32630)

cat("-> Re-projecting and precision-clipping raster to eliminate edge offsets...\n")
# 1. Project the raw cassava map to UTM first to establish the clean target grid
cassava_utm_full <- terra::project(cassavaMap, "EPSG:32630")
# 2. Crop and Mask AFTER reprojection so the raster pixels match vector cuts exactly
cassava_bg_utm <- terra::crop(cassava_utm_full, country_borders_utm) %>% 
  terra::mask(country_borders_utm)

#-----------------------------------------------------------------------------
# 2. Map Field-Confirmed Points to the Host Landscape ----
#-----------------------------------------------------------------------------
cat(">>> Mapping field-confirmed epicenters to landscape nodes...\n")

true_positives <- dataEpi %>%
  filter(`EACMV-Ug` == "yes") %>%
  mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
  filter(!is.na(Longitude), !is.na(Latitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(32630) 

nodes_sf             <- st_as_sf(empirical_nodes, coords = c("x", "y"), crs = 32630)
nearest_node_indices <- st_nearest_feature(true_positives, nodes_sf) %>% unique()

#-----------------------------------------------------------------------------
# 3. Execute Unified Connected Network Forward Simulation ----
#-----------------------------------------------------------------------------
cat(">>> Running aligned network epidemic forward projection over host matrix...\n")

simulation_history <- vector("list", total_years)
status             <- rep("Susceptible", N_nodes)
status[nearest_node_indices] <- "Infected"

simulation_history[[1]] <- empirical_nodes %>%
  mutate(Year = 1, Status = status)

for (year in 2:total_years) {
  current_infected    <- which(status == "Infected")
  current_susceptible <- which(status == "Susceptible")
  
  if (length(current_infected) > 0 && length(current_susceptible) > 0) {
    # Extract the slice connecting susceptible individuals to active infective nodes
    sub_matrix <- env_matrix[current_susceptible, current_infected, drop = FALSE]
    
    # Vectorized connection-link joint probability calculation
    escape_probabilities <- exp(rowSums(log(1 - pmin(sub_matrix * transmission_scale_val, 0.9999))))
    infection_risk       <- 1 - escape_probabilities
    
    # Stochastic transition assignment
    new_infections <- current_susceptible[runif(length(current_susceptible)) < infection_risk]
    
    if (length(new_infections) > 0) {
      status[new_infections] <- "Infected"
    }
  }
  
  simulation_history[[year]] <- empirical_nodes %>%
    mutate(Year = year, Status = status)
}

animation_df <- bind_rows(simulation_history) %>%
  mutate(
    risk_norm     = (risk - min(risk)) / (max(risk) - min(risk)), 
    display_value = if_else(Status == "Susceptible", risk_norm, risk_norm + 2)
  )

infected_only_df <- animation_df %>% filter(Status == "Infected")

#-----------------------------------------------------------------------------
# 4. Render Spatial Spatiotemporal Progression Animation ----
#-----------------------------------------------------------------------------
cat(">>> Compiling map visualization steps via gganimate...\n")

anim_plot <- ggplot() +
  # LAYER 1: Corrected continuous Green-to-White base raster
  geom_spatraster(data = cassava_bg_utm, maxcell = 150000, alpha = 0.85) +
  scale_fill_gradient(
    low = "#f7fcf5",   
    high = "#00441b",  
    na.value = "transparent",
    name = "Cassava Production Density",
    labels = scales::label_comma()
  ) +
  
  ggnewscale::new_scale_color() +
  
  # LAYER 2: Infection outbreak markers (Significantly reduced dot sizes)
  geom_point(
    data = infected_only_df,
    aes(x = x, y = y, color = display_value, size = risk), 
    alpha = 0.95
  ) +
  scale_color_gradientn(
    colors = c("#fee0d2", "#fc9272", "#ef3b2c", "#cb181d"), 
    values = c(0.0, 0.35, 0.70, 1.0),
    breaks = c(2.2, 2.8),
    labels = c("Initial Entry", "Severe Outbreak"),
    name = "Invasion Status"
  ) +
  scale_size_continuous(range = c(0.1, 1.0), name = "Outbreak Risk Magnitude") +
  
  # LAYER 3: National Border Overlays - Adjusted to be clean and light
  geom_sf(
    data = country_borders_utm, 
    fill = NA, 
    color = "#95a5a6",  # Light slate grey
    linewidth = 0.25,   # Thinned out for crisp alignment limits
    alpha = 0.60, 
    inherit.aes = FALSE
  ) +
  
  annotation_north_arrow(
    location = "tl", which_north = "true",
    pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
    style = north_arrow_fancy_orienteering(fill = c("grey20", "white"))
  ) +
  annotation_scale(
    location = "bl", width_hint = 0.25,
    unit_category = "metric", style = "bar" 
  ) +
  
  theme_minimal(base_family = "sans") + 
  labs(
    title = "West Africa: 30-Year Calibrated Network SI Progression Simulation",
    subtitle = "Year: {frame_time} of 30 | Base = Production Continuum | Borders = National Limits",
    x = "Easting (UTM 30N / Meters)",
    y = "Northing"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "#1a1a1a"),
    plot.subtitle = element_text(size = 12, color = "grey20", margin = margin(b = 10)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
    panel.background = element_rect(fill = "#fafafa", color = NA)
  ) +
  transition_time(Year) +
  ease_aes('linear')

#-----------------------------------------------------------------------------
# High-Quality Video Rendering Output Configuration ----
#-----------------------------------------------------------------------------
cat("-> Compiling video frames...\n")

rendered_anim <- gganimate::animate(
  plot     = anim_plot,
  nframes  = total_years,
  fps      = 1,                                   
  res      = 150,                                
  width    = 1200,                              
  height   = 850,
  renderer = gifski_renderer() 
)

anim_save(
  filename  = "outputs/grounded_epidemic_progression_highres.gif",
  animation = rendered_anim
)

gganimate::animate(
  plot      = anim_plot,
  nframes   = total_years,
  fps       = 1,
  res       = 150,
  width     = 1200,
  height    = 850,
  renderer  = av_renderer(file = "outputs/grounded_epidemic_progression_highres.mp4")
)

cat("✔ Done! Post-processed with crisp network alignments and optimized visual sizes.\n")

# #=============================================================================
# # scripts/simulate_progression.R
# # Grounded 30-Year Spatiotemporal Outbreak Projection & Video Generation
# #=============================================================================
# 
# library(tidyverse)
# library(terra)
# library(sf)
# library(ggplot2)
# library(tidyterra)
# library(gganimate)
# library(ggspatial) 
# library(av)
# 
# #-----------------------------------------------------------------------------
# # 1. Ingest Pipeline Spatial Assets & Parameters ----
# #-----------------------------------------------------------------------------
# if(!file.exists("outputs/nodes.rds") || !file.exists("outputs/distance_matrix.rds")) {
#   stop("Missing structural inputs. Please run 'scripts/network_nodes.R' first.")
# }
# if(!file.exists("outputs/optimal_abc_parameters.rds")) {
#   stop("Parameters uncalibrated. Please run 'scripts/ibm_smc_abc.R' first.")
# }
# 
# # Load the empirical host landscape matrix generated by network_nodes.R
# empirical_nodes <- readRDS("outputs/nodes.rds")
# empirical_dist_matrix <- readRDS("outputs/distance_matrix.rds")
# 
# # Load ABC optimized parameter metrics
# abc_params <- readRDS("outputs/optimal_abc_parameters.rds")
# beta_val  <- abc_params %>% filter(parameter == "beta") %>% pull(optimal_value)
# gamma_val <- abc_params %>% filter(parameter == "gamma_dist") %>% pull(optimal_value)
# 
# N_nodes <- nrow(empirical_nodes)
# total_years <- 30
# 
# # ============================================================================
# # TARGET FIX: Sourcing & Realignment for Perfect Border Edge Snapping
# # ============================================================================
# cat(">>> Ingesting pre-cropped landscape boundaries from cassMap.R...\n")
# source("scripts/cassMap.R") 
# 
# cat("-> Re-projecting vector assets to UTM 30N (EPSG:32630)...\n")
# country_borders_utm <- st_transform(three_countries_proj, 32630)
# 
# cat("-> Re-projecting and precision-clipping raster to eliminate edge offsets...\n")
# # 1. Project the raw cassava map to UTM first to establish the clean target grid
# cassava_utm_full <- terra::project(cassavaMap, "EPSG:32630")
# # 2. Crop and Mask AFTER reprojection so the raster pixels match vector cuts exactly
# cassava_bg_utm <- terra::crop(cassava_utm_full, country_borders_utm) %>% 
#   terra::mask(country_borders_utm)
# 
# #-----------------------------------------------------------------------------
# # 2. Map Field-Confirmed Points to the Host Landscape ----
# #-----------------------------------------------------------------------------
# cat(">>> Mapping field-confirmed epicenters to landscape nodes...\n")
# 
# true_positives <- dataEpi %>%
#   filter(`EACMV-Ug` == "yes") %>%
#   mutate(Longitude = as.numeric(Longitude), Latitude = as.numeric(Latitude)) %>%
#   filter(!is.na(Longitude), !is.na(Latitude)) %>%
#   st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
#   st_transform(32630) 
# 
# nodes_sf <- st_as_sf(empirical_nodes, coords = c("x", "y"), crs = 32630)
# nearest_node_indices <- st_nearest_feature(true_positives, nodes_sf) %>% unique()
# 
# #-----------------------------------------------------------------------------
# # 3. Execute Grounded Forward Projection Simulation ----
# #-----------------------------------------------------------------------------
# cat(">>> Running forward epidemic projection over host matrix...\n")
# 
# simulation_history <- vector("list", total_years)
# status <- rep("Susceptible", N_nodes)
# status[nearest_node_indices] <- "Infected"
# 
# simulation_history[[1]] <- empirical_nodes %>%
#   mutate(Year = 1, Status = status)
# 
# for (year in 2:total_years) {
#   current_infected <- which(status == "Infected")
#   current_susceptible <- which(status == "Susceptible")
#   
#   if (length(current_infected) > 0 && length(current_susceptible) > 0) {
#     escape_probabilities <- sapply(current_susceptible, function(s_idx) {
#       distances <- empirical_dist_matrix[s_idx, current_infected]
#       kernel = beta_val / (1 + (distances / 1000)^gamma_val)
#       prod(1 - kernel)
#     })
#     
#     infection_risk <- 1 - escape_probabilities
#     new_infections <- current_susceptible[runif(length(current_susceptible)) < infection_risk]
#     
#     if (length(new_infections) > 0) {
#       status[new_infections] <- "Infected"
#     }
#   }
#   
#   simulation_history[[year]] <- empirical_nodes %>%
#     mutate(Year = year, Status = status)
# }
# 
# animation_df <- bind_rows(simulation_history) %>%
#   mutate(
#     risk_norm = (risk - min(risk)) / (max(risk) - min(risk)), 
#     display_value = if_else(Status == "Susceptible", risk_norm, risk_norm + 2)
#   )
# 
# infected_only_df <- animation_df %>% filter(Status == "Infected")
# 
# #-----------------------------------------------------------------------------
# # 4. Render Spatial Spatiotemporal Progression Animation ----
# #-----------------------------------------------------------------------------
# cat(">>> Compiling map visualization steps via gganimate...\n")
# 
# anim_plot <- ggplot() +
#   # LAYER 1: Corrected continuous Green-to-White base raster
#   geom_spatraster(data = cassava_bg_utm, maxcell = 150000, alpha = 0.85) +
#   scale_fill_gradient(
#     low = "#f7fcf5",   
#     high = "#00441b",  
#     na.value = "transparent",
#     name = "Cassava Production Density",
#     labels = scales::label_comma()
#   ) +
#   
#   ggnewscale::new_scale_color() +
#   
#   # LAYER 2: Infection outbreak markers (Significantly reduced dot sizes)
#   geom_point(
#     data = infected_only_df,
#     aes(x = x, y = y, color = display_value, size = risk), 
#     alpha = 0.95
#   ) +
#   scale_color_gradientn(
#     colors = c("#fee0d2", "#fc9272", "#ef3b2c", "#cb181d"), 
#     values = c(0.0, 0.35, 0.70, 1.0),
#     breaks = c(2.2, 2.8),
#     labels = c("Initial Entry", "Severe Outbreak"),
#     name = "Invasion Status"
#   ) +
#   # ENHANCEMENT: Dots shrunk from c(0.4, 2.2) to c(0.1, 1.0) for a precise tracking overview
#   scale_size_continuous(range = c(0.1, 1.0), name = "Outbreak Risk Magnitude") +
#   
#   # LAYER 3: National Border Overlays - Adjusted to be clean and light
#   geom_sf(
#     data = country_borders_utm, 
#     fill = NA, 
#     color = "#95a5a6",  # Lighter slate grey
#     linewidth = 0.25,   # Thinned out from 0.35 to 0.25 for crisp alignment limits
#     alpha = 0.60, 
#     inherit.aes = FALSE
#   ) +
#   
#   annotation_north_arrow(
#     location = "tl", which_north = "true",
#     pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
#     style = north_arrow_fancy_orienteering(fill = c("grey20", "white"))
#   ) +
#   annotation_scale(
#     location = "bl", width_hint = 0.25,
#     unit_category = "metric", style = "bar" 
#   ) +
#   
#   theme_minimal(base_family = "sans") + 
#   labs(
#     title = "West Africa: 30-Year Grounded CMD/UG Invasion Simulation",
#     subtitle = "Year: {frame_time} of 30 | Base = Production Continuum | Borders = National Limits",
#     x = "Easting (UTM 30N / Meters)",
#     y = "Northing"
#   ) +
#   theme(
#     plot.title = element_text(face = "bold", size = 16, color = "#1a1a1a"),
#     plot.subtitle = element_text(size = 12, color = "grey20", margin = margin(b = 10)),
#     legend.position = "right",
#     legend.title = element_text(face = "bold", size = 9),
#     panel.grid.major = element_line(color = "grey92", size = 0.25),
#     panel.background = element_rect(fill = "#fafafa", color = NA)
#   ) +
#   transition_time(Year) +
#   ease_aes('linear')
# 
# #-----------------------------------------------------------------------------
# # High-Quality Video Rendering Output Configuration ----
# #-----------------------------------------------------------------------------
# cat("-> Compiling video frames...\n")
# 
# rendered_anim <- gganimate::animate(
#   plot     = anim_plot,
#   nframes  = total_years,
#   fps      = 1,                                 
#   res      = 150,                               
#   width    = 1200,                              
#   height   = 850,
#   renderer = gifski_renderer() 
# )
# 
# anim_save(
#   filename  = "outputs/grounded_epidemic_progression_highres.gif",
#   animation = rendered_anim
# )
# 
# gganimate::animate(
#   plot      = anim_plot,
#   nframes   = total_years,
#   fps       = 1,
#   res       = 150,
#   width     = 1200,
#   height    = 850,
#   renderer  = av_renderer(file = "outputs/grounded_epidemic_progression_highres.mp4")
# )
# 
# cat("✔ Done! Post-processed with crisp alignments and optimized visual sizes.\n")