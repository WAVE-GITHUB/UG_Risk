#======================================================
# Data wraggling and Distance Matrix computation
#======================================================

source("R/library.R")

dataEpi <- read_xlsx("data/epidata.xlsx", sheet = "Field & Lab")

library(readxl)
library(dplyr)
library(sf)

# 1. Load and Clean Epidemiological Data ----
epi_clean <-dataEpi %>%
  # Convert types and filter out invalid coordinates early
  mutate(
    Latitude  = as.numeric(Latitude),
    Longitude = as.numeric(Longitude)
  ) %>%
  filter(!is.na(Latitude), !is.na(Longitude)) %>% 
  
  # Recode diagnostic variables
  mutate(
    CMD   = if_else(CMD_Severity > 1, 1, 0),
    UG    = if_else(`EACMV-Ug` == "yes", 1, 0),
    Virus = if_else(Diagnostic != "Healthy", 1, 0)
  ) %>% 
  
  # Aggregate by field
  group_by(`Full _Field_ID`) %>%
  summarise(
    # Spatial anchors
    lon = mean(Longitude),
    lat = mean(Latitude),
    
    # Sample size & epidemiological metrics
    nPlants  = n(),
    CMD_prev = mean(CMD, na.rm = TRUE),
    UG_prev  = mean(UG, na.rm = TRUE),
    wf_mean  = mean(Total_Whitefly_Count, na.rm = TRUE),
    
    # Contextual metadata
    Country     = first(Country),
    Admin_1     = first(Admin_1),
    Admin_2     = first(Admin_2),
    Survey_Year = first(Survey_Year),
    .groups = "drop"
  ) %>% 
  
  # Determine field-level infection status
  mutate(infected = if_else(UG_prev > 0, 1, 0)) %>% 
  
  # Convert to simple features (WGS84) and project (UTM Zone 30N)
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  st_transform(crs = 32630) %>% 
  
  # Append projected coordinates back into the dataframe cleanly
  mutate(
    x = st_coordinates(.)[, 1],
    y = st_coordinates(.)[, 2]
  )

# 2. Compute Inter-Farm Distance Matrix (in km) ----
# Drops the 'm' or 'km' units class to ensure raw numeric matrix compatibility
dist_mat <- as.matrix(st_distance(epi_clean)) / 1000
units(dist_mat) <- NULL
