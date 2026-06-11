#===========================================================
# riskMap.R
# Gravity-Weighted Landscape Invasion Risk Surface
#===========================================================

library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(tidyterra)

#-----------------------------------------------------------
# 1. Environment & Path Prerequisites
#-----------------------------------------------------------
# Ensure directories exist
if(!dir.exists("outputs")) dir.create("outputs")

# [CRITICAL] Assumes the following layers are already loaded from your pipeline:
# - dataEpi (The foundational excel data frame)
# - cassava_3countries (The cropped/masked SpatRaster)
# - road_density_norm (The normalized focal road density SpatRaster)
# - major_roads (The original sf line geometry dataframe)

#-----------------------------------------------------------
# 2. Project Base Grids to Metric Coordinates (UTM 30N)
#-----------------------------------------------------------
# Forcing analysis into meters ensures distance kernels evaluate accurately
target_crs <- st_crs(32630)$wkt

if(crs(cassava_3countries, proj = TRUE) != crs(target_crs, proj = TRUE)) {
  cassava_utm <- project(cassava_3countries, target_crs, method = "bilinear")
} else {
  cassava_utm <- cassava_3countries
}

# Sync the road grid to match the metric raster template exactly
road_risk_utm <- project(road_density_norm, cassava_utm, method = "bilinear")

#-----------------------------------------------------------
# 3. Extract and Project Pathogen Positive Tracking Nodes
#-----------------------------------------------------------
ug_pos <- dataEpi %>%
  filter(`EACMV-Ug` == "yes") %>%
  mutate(
    Longitude = as.numeric(Longitude), 
    Latitude = as.numeric(Latitude)
  ) %>% 
  filter(!is.na(Longitude), !is.na(Latitude))

ug_sf <- st_as_sf(ug_pos, coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(32630)

infected_farms_vect <- vect(ug_sf)

#-----------------------------------------------------------
# 4. Compute Regional Long-Distance Dispersal Kernel
#-----------------------------------------------------------
# Calculate continuous distance from every cell to nearest outbreak node (in meters)
distance_to_infected <- distance(cassava_utm, infected_farms_vect)

# Parameters: Using 80km sigma to realistically capture regional transit corridors
alpha <- 2.0
sigma_trade <- 80000  

kernel_raster <- 1 / (1 + (distance_to_infected / sigma_trade)^alpha)
kernel_raster_norm <- kernel_raster / global(kernel_raster, "max", na.rm = TRUE)[[1]]

#-----------------------------------------------------------
# 5. Build the Non-Linear Multiplicative Gravity Index
#-----------------------------------------------------------
# FIX: Exponentiating the layers amplifies contrast between high and low values.
# This prevents high-production hotspots from getting flattened by the math.
host_weight         <- 1.5  # Amplifies the penalty for dense host presence
connectivity_weight <- 1.2  # Weights major transit routes

weighted_base_risk <- (cassava_utm^host_weight) * (road_risk_utm^connectivity_weight)

# Normalize structural layer (0 to 1)
structural_risk <- weighted_base_risk / global(weighted_base_risk, "max", na.rm = TRUE)[[1]]

# Combine: Outbreak Proximity x (Host Density x Trade Connectivity)
final_weighted_risk <- kernel_raster_norm * structural_risk

# Strip out structural zeros to isolate active agricultural zones
final_weighted_risk <- mask(final_weighted_risk, cassava_utm)
final_weighted_risk[final_weighted_risk <= 0] <- NA
names(final_weighted_risk) <- "Weighted_Invasion_Risk"

#-----------------------------------------------------------
# 6. Filter Regional Road Vectors for Clean Map Overlay
#-----------------------------------------------------------
if(exists("major_roads")) {
  display_roads <- major_roads %>%
    filter(road_type %in% c("motorway", "trunk", "primary", "secondary")) %>%
    st_transform(32630)
}

#-----------------------------------------------------------
# 7. High-Contrast Render & Export Engine
#-----------------------------------------------------------
risk_map_plot <- ggplot() +
  # Draw risk surface with a vivid yellow-to-red gradient scale
  geom_spatraster(data = final_weighted_risk, aes(fill = Weighted_Invasion_Risk)) +
  
  scale_fill_gradientn(
    name = "Risk Intensity",
    colors = c("#ffeb3b", "#ff9800", "#f44336", "#b71c1c"), # Clear heat progression
    na.value = "transparent",
    trans = "sqrt",  # Square-root scale balances variance without crushing high-production peaks
    labels = scales::label_scientific(digits = 2)
  ) +
  
  # Overlay primary transportation networks
  geom_sf(data = display_roads, color = "black", alpha = 0.20, size = 0.2) +
  
  # Overlay known infected epicenters
  geom_sf(data = st_as_sf(ug_sf), color = "white", fill = "#00ffff", shape = 21, size = 2.0, stroke = 0.7) +
  
  # Layout thematic adjustments
  theme_minimal() +
  labs(
    title = "West Africa: CMD/UG Spatial Invasion Risk Surface",
    subtitle = "Gravity-Weighted Host Clusters & Trade Connectivity Corridors",
    x = "Easting (UTM 30N)",
    y = "Northing"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8),
    panel.grid.major = element_line(color = "grey92", size = 0.2),
    plot.background = element_rect(fill = "white", color = NA)
  )

# Display to the active screen device
print(risk_map_plot)

# Export high-resolution file for publication
ggsave(
  filename = "outputs/UG_HighContrast_RiskMap.png",
  plot = risk_map_plot,
  width = 11,
  height = 7,
  dpi = 300
)

risk_map <- final_weighted_risk