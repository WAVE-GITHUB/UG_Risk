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

#-----------------------------------------------------------
# Merge
#-----------------------------------------------------------

roads_all <- rbind(
  roads_civ,
  roads_sle,
  roads_gin
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
