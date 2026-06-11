#===========================================================
# Cassava Map
#===========================================================

# ===========================================================
# 1. Load Cassava Map (Keep the raster object intact!)
# ===========================================================
# terra::rast() replaces raster::raster()
cassavaMap <- rast("data/CassavaMap_Prod_v1.tif") 

# ===========================================================
# 2. Load and Combine Shapefiles
# ===========================================================
civ     <- st_read("data/countriesMaps/civ_admin/civ_admin1_em.shp")
sierraL <- st_read("data/countriesMaps/sle_admin/sle_admin1.shp")
guinea  <- st_read("data/countriesMaps/gna_admin/gis_osm_adminareas_a_free_1.shp")

civ2 <- civ %>% 
  dplyr::select(country = adm0_name, admin1 = adm1_name)

sle2 <- sierraL %>% 
  dplyr::select(country = adm0_name, admin1 = adm1_name)

gin2 <- guinea %>% 
  transmute(country = "Guinea", admin1 = name)

# Combine, fix geometries, and dissolve boundaries
three_countries <- rbind(civ2, sle2, gin2) %>% 
  st_make_valid() %>%
  st_union()

# ===========================================================
# 3. Project Vector Data to Match Raster CRS
# ===========================================================
# Project the sf object directly using the raster's CRS
three_countries_proj <- st_transform(three_countries, crs(cassavaMap))

# Convert the sfc geometry back into a proper sf dataframe
three_countries_sf <- st_as_sf(three_countries_proj)

# ===========================================================
# 4. Crop and Mask (Native terra workflow)
# ===========================================================
# terra::crop accepts sf objects directly without converting to 'Spatial'
cassava_crop <- crop(cassavaMap, three_countries_proj)

# Mask the cropped raster to the exact vector boundaries
cassava_3countries <- mask(cassava_crop, three_countries_sf)