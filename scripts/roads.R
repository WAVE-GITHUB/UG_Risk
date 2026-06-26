#===========================================================
# roads.R
# Road accessibility layer
#===========================================================


#-----------------------------------------------------------
# Read roads
#-----------------------------------------------------------

roads_civ <- st_read(
  "data/roadMaps/civ_roads/roads_lines.shp",
  quiet = TRUE
)

roads_sle <- st_read(
  "data/roadMaps/sle_roads/roads_lines.shp",
  quiet = TRUE
)

roads_gin <- st_read(
  "data/roadMaps/gna_roads/roads_lines.shp",
  quiet = TRUE
)

roads_ben <- st_read(
  "data/roadMaps/ben_roads/roads_lines.shp",
                     quiet = TRUE)

roads_bfa <- st_read(
  "data/roadMaps/bfa_roads/roads_lines.shp",
                     quiet = TRUE)

roads_gha <- st_read(
  "data/roadMaps/gha_roads/roads_lines.shp",
                     quiet = TRUE)

roads_lbr <- st_read(
  "data/roadMaps/lbr_roads/roads_lines.shp",
                     quiet = TRUE)

roads_ner <- st_read(
  "data/roadMaps/ner_roads/roads_lines.shp",
                     quiet = TRUE)

roads_nga <- st_read(
  "data/roadMaps/nga_roads/hotosm_nga_roads_lines_shp.shp",
                     quiet = TRUE)

roads_tgo <- st_read(
  "data/roadMaps/tgo_roads/roads_lines.shp",
                     quiet = TRUE)

#-----------------------------------------------------------
# Harmonize schema
#-----------------------------------------------------------

roads_civ <- roads_civ %>%
  transmute(
    country = "Cote d'Ivoire",
    road_type = highway,
    geometry
  )

roads_sle <- roads_sle %>%
  transmute(
    country = "Sierra Leone",
    road_type = highway,
    geometry
  )

roads_gin <- roads_gin %>%
  transmute(
    country = "Guinea",
    road_type = highway,
    geometry
  )

roads_ben <-  roads_ben  %>% 
  transmute(
    country = "Benin",
    road_type = highway,
    geometry
  )

roads_bfa <-  roads_bfa %>% 
  transmute(
    country = "Burkina-Faso",
    road_type = highway,
    geometry
  )

roads_gha <- roads_gha %>% 
  transmute(
    country = "Ghana",
    road_type = highway,
    geometry
  )

roads_lbr <-  roads_lbr %>% 
  transmute(
    country = "Liberia",
    road_type = highway,
    geometry
  )

roads_ner <-  roads_ner %>% 
  transmute(
    country = "Niger",
    road_type = highway,
    geometry
  )

roads_nga <-  roads_nga %>% 
  transmute(
    country = "Nigeria",
    road_type = highway,
    geometry
  )

roads_tgo <-  roads_tgo %>% 
  transmute(
    country = "Togo",
    road_type = highway,
    geometry
  )
#-----------------------------------------------------------
# Merge
#-----------------------------------------------------------

roads_all <- rbind(
  roads_civ,
  roads_sle,
  roads_gin, 
  roads_ben, 
  roads_bfa,
  roads_gha,
  roads_lbr,
  roads_ner,
  roads_nga,
  roads_tgo
)

#-----------------------------------------------------------
# Keep epidemiologically relevant roads
#-----------------------------------------------------------

major_roads <- roads_all %>%
  filter(
    road_type %in% c(
      "motorway",
      "trunk",
      "primary",
      "secondary",
      "tertiary",
      "unclassified",
      "residential"
    )
  )

#-----------------------------------------------------------
# Remove invalid geometries
#-----------------------------------------------------------

major_roads <- major_roads %>%
  st_make_valid()

#-----------------------------------------------------------
# Reproject to cassava raster CRS
#-----------------------------------------------------------

major_roads <- st_transform(
  major_roads,
  crs(cassava_3countries)
)

#-----------------------------------------------------------
# Road density raster
#-----------------------------------------------------------

road_vect <- vect(major_roads)

road_density <- rasterize(
  road_vect,
  cassava_3countries,
  field = 1,
  fun = "sum",
  background = 0
)

#-----------------------------------------------------------
# Smooth density
#-----------------------------------------------------------

road_density <- focal(
  road_density,
  w = 5,
  fun = mean,
  na.policy = "omit"
)

#-----------------------------------------------------------
# Scale 0-1
#-----------------------------------------------------------

mn <- global(road_density, "min", na.rm = TRUE)[[1]]
mx <- global(road_density, "max", na.rm = TRUE)[[1]]

road_density_norm <- (road_density - mn) / (mx - mn)

road_risk <- road_density_norm
