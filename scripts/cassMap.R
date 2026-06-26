#===========================================================
# Cassava Map - Country Borders (Admin 0)
#===========================================================

# ===========================================================
# 1. Load Cassava Map (Keep the raster object intact!)
# ===========================================================
cassavaMap <- rast("data/CassavaMap_Prod_v1.tif") 

# ===========================================================
# 2. Load and Combine Shapefiles
# ===========================================================
civ     <- st_read("data/countriesMaps/civ_admin/civ_admin1_em.shp")
sierraL <- st_read("data/countriesMaps/sle_admin/sle_admin1.shp")
guinea  <- st_read("data/countriesMaps/gna_admin/gis_osm_adminareas_a_free_1.shp")
ben <- st_read("data/countriesMaps/ben_admin/geoBoundaries-BEN-ADM3_simplified.shp")
bfa <- st_read("data/countriesMaps/bfa_admin/geoBoundaries-BFA-ADM2_simplified.shp")
gha <- st_read("data/countriesMaps/gha_admin/geoBoundaries-GHA-ADM2_simplified.shp")
lbr <-st_read("data/countriesMaps/lbr_admin/lbr_admin0_em.shp")
ner <-st_read("data/countriesMaps/ner_admin/geoBoundaries-NER-ADM2.shp")
nga <- st_read("data/countriesMaps/nga_admin/geoBoundaries-NGA-ADM2_simplified.shp")
tgo <- st_read("data/countriesMaps/tgo_admin/geoBoundaries-TGO-ADM2_simplified.shp")


civ2 <- civ %>% 
  dplyr::select(country = adm0_name, admin1 = adm1_name)

sle2 <- sierraL %>% 
  dplyr::select(country = adm0_name, admin1 = adm1_name)

gin2 <- guinea %>% 
  transmute(country = "Guinea", admin1 = name)

ben2 <- ben %>% 
  transmute(country = "Benin", admin1 = shapeName)

bfa2 <- bfa %>% 
  transmute(country = "Burkina-Faso", admin1 = shapeName)

gha2 <- gha %>% 
  transmute(country = "Ghana", admin1 = shapeName)

lbr2 <- lbr %>% 
  dplyr::select(country = adm0_name,  admin1 = adm0_name1)

ner2 <- ner %>% 
  transmute(country = "Niger", admin1 = shapeName)

nga2 <- nga %>% 
  transmute(country = "Nigeria", admin1 = shapeName)

tgo2 <- tgo %>% 
  transmute(country = "Togo", admin1 = shapeName)

# FIX: Group by country before summarizing to keep Admin 0 national borders intact
three_countries <- rbind(civ2, sle2, gin2, ben2, bfa2, gha2, lbr2, ner2, nga2,tgo2) %>% 
  st_make_valid() %>%
  group_by(country) %>% 
  summarize(.groups = "drop") # This dissolves Admin 1 into distinct Admin 0 borders

# ===========================================================
# 3. Project Vector Data to Match Raster CRS
# ===========================================================
three_countries_proj <- st_transform(three_countries, crs(cassavaMap))

# ===========================================================
# 4. Crop and Mask (Native terra workflow)
# ===========================================================
cassava_crop       <- crop(cassavaMap, three_countries_proj)
cassava_3countries <- mask(cassava_crop, three_countries_proj)