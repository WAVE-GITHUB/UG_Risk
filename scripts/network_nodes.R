############################################################
## Build epidemiological network from cassava suitability
## and road-mediated movement (Border-Weighted Network)
############################################################


# Ensure output repository structure is active
if(!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

# CONFIGURATION PARAMETER: Cross-border transmission weight modifier
# 1.0 = No penalty (Original gravity model)
# 0.1 = Cross-border movement is penalized by 90% (forcing high domestic density)
border_weight_modifier <- 0.15 

#-----------------------------------------------------------
# 1. Select and Aggregate Occupied Cassava Cells ----
#-----------------------------------------------------------
cat("Aggregating risk grid to regional landscape nodes...\n")

risk_agg <- aggregate(
  risk_map,
  fact = 10,
  fun = mean,
  na.rm = TRUE
)

# Extract non-NA geometric points from the aggregated spatial model
nodes_spatial <- as.points(risk_agg)

# Convert to an sf object to perform a spatial join with country boundaries
nodes_sf <- st_as_sf(nodes_spatial)
st_crs(nodes_sf) <- 32630

# Standardize risk column naming
colnames(nodes_sf)[1] <- "risk"

# Clean baseline records
nodes_sf <- nodes_sf %>%
  filter(!is.na(risk), risk > 0)

#-----------------------------------------------------------
# NEW STEP: Spatial Intersection with Country Boundaries ----
#-----------------------------------------------------------
cat("Intersecting network nodes with country perimeters...\n")

# Verify study_borders exists from your previous geometry processing step
if(!exists("study_borders")) {
  stop("CRITICAL DEPENDENCY MISSING: 'study_borders' spatial object was not found.")
}

# Perform spatial point-in-polygon join to tag each node with its respective country name
nodes_mapped_sf <- st_join(nodes_sf, study_borders, join = st_intersects)

# Extract coordinates and convert back to a clean modeling data frame
nodes_df <- nodes_mapped_sf %>%
  mutate(
    node_id = row_number(),
    x = st_coordinates(.)[,1],
    y = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry() %>%
  # Fallback allocation in case a border-edge node sits pixel-wise outside the vector outlines
  mutate(country = ifelse(is.na(country), "Border-Edge-Artifact", country)) %>%
  dplyr::select(node_id, x, y, risk, country)

cat("Successfully generated network nodes. Count:", nrow(nodes_df), "\n")
print(table(nodes_df$country)) # Prints node count per country for diagnostics

saveRDS(nodes_df, "outputs/nodes.rds")

#-----------------------------------------------------------
# 2. Distance Matrix Calculation ----
#-----------------------------------------------------------
cat("Computing geographical Euclidean distance matrix...\n")

coords_matrix <- as.matrix(nodes_df[, c("x", "y")])
distance_matrix <- as.matrix(dist(coords_matrix))

rownames(distance_matrix) <- nodes_df$node_id
colnames(distance_matrix) <- nodes_df$node_id

saveRDS(distance_matrix, "outputs/distance_matrix.rds")

#-----------------------------------------------------------
# 3. Compute Gravity-Based Movement with Border Penalties ----
#-----------------------------------------------------------
cat("Formulating border-weighted gravity network connectivity matrix...\n")

risk_vector <- nodes_df$risk
risk_product_matrix <- outer(risk_vector, risk_vector, "*")

# Base landscape gravity structure layer
base_movement_matrix <- risk_product_matrix / (1 + (distance_matrix)^2)

# NEW MATH: Generate Country Match Multiplier Matrix
cat("Applying international border restrictions and friction modifiers...\n")
country_vector <- nodes_df$country

# Matrix outer match: Returns TRUE if Row Node and Column Node share the same country string
country_match_matrix <- outer(country_vector, country_vector, "==")

# Construct modifier array: 1.0 for matches, border_weight_modifier for international links
border_penalty_matrix <- ifelse(country_match_matrix, 1.0, border_weight_modifier)

# Apply the border weight drop to the base network weights
movement_matrix <- base_movement_matrix * border_penalty_matrix

# Clear diagonal tracking loops
diag(movement_matrix) <- 0

saveRDS(movement_matrix, "outputs/movement_matrix.rds")
cat("✔ Border-weighted network execution complete. Matrices saved to /outputs/.\n")