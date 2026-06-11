############################################################
## Build epidemiological network from cassava suitability
## and road-mediated movement
############################################################

library(terra)
library(sf)
library(tidyverse)
library(igraph)

# Ensure output repository structure is active
if(!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

#-----------------------------------------------------------
# 1. Select and Aggregate Occupied Cassava Cells ----
#-----------------------------------------------------------
# High-resolution cell extraction can lead to millions of network nodes,
# causing exponential memory exhaustion (O(N^2)) during matrix generation.
# We aggregate by factor = 10 to establish clean regional landscape patches.

cat("Aggregating risk grid to regional landscape nodes...\n")

risk_agg <- aggregate(
  risk_map,
  fact = 10,
  fun = mean,
  na.rm = TRUE
)

# Extract non-NA geometric points from the aggregated spatial model
nodes_spatial <- as.points(risk_agg)

# Convert cleanly to a data frame containing explicit x/y coordinates
nodes_df <- as.data.frame(nodes_spatial, geom = "XY")

# Standardize names safely regardless of what the input raster layer was named
colnames(nodes_df)[1] <- "risk"

# Clean out absolute zeros, negative noise, and missing value records
nodes_df <- nodes_df %>%
  filter(!is.na(risk), risk > 0) %>%
  mutate(node_id = row_number()) %>%
  select(node_id, x = x, y = y, risk)

cat("Successfully generated network nodes. Count:", nrow(nodes_df), "\n")

# Save primary nodes database tracking file
saveRDS(
  nodes_df,
  "outputs/nodes.rds"
)

#-----------------------------------------------------------
# 2. Distance Matrix Calculation ----
#-----------------------------------------------------------
cat("Computing geographical Euclidean distance matrix...\n")

# Matrix operations require a clean coordinate pair matrix
coords_matrix <- as.matrix(nodes_df[, c("x", "y")])

# Calculate distances between all node combinations
distance_matrix <- as.matrix(dist(coords_matrix))

# Apply matching node IDs to matrix row and column headers for network tracking
rownames(distance_matrix) <- nodes_df$node_id
colnames(distance_matrix) <- nodes_df$node_id

saveRDS(
  distance_matrix,
  "outputs/distance_matrix.rds"
)

#-----------------------------------------------------------
# 3. Compute Gravity-Based Movement / Adjacency Matrix ----
#-----------------------------------------------------------
cat("Formulating gravity-based network connectivity matrix...\n")

# Constructing the actual movement connectivity weight matrices
# Formula: Movement(i,j) = (Risk_i * Risk_j) / (1 + Distance_ij^2)
risk_vector <- nodes_df$risk

# Outer product multiplies vector values across pairs: Risk_i * Risk_j
risk_product_matrix <- outer(risk_vector, risk_vector, "*")

# Apply an inverse-distance power decay function to capture road-mediated movement
# Adding 1 prevents division by zero anomalies on the diagonal axis
movement_matrix <- risk_product_matrix / (1 + (distance_matrix)^2)

# Set self-loop interactions to absolute zero
diag(movement_matrix) <- 0

saveRDS(
  movement_matrix,
  "outputs/movement_matrix.rds"
)

cat("All epidemiological network matrices saved successfully to /outputs/.\n")

