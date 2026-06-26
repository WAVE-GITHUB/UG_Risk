#=============================================================================
# scripts/compute_economic_losses.R
# Bio-Economic Valuation & Spatial Damage Function Processor
#=============================================================================

library(terra)
library(dplyr)
library(tidyr)

#-----------------------------------------------------------------------------
# 1. Load the Harvest Map Asset & Price Parameters
#-----------------------------------------------------------------------------
cat(">>> Ingesting base crop harvest raster metrics...\n")

# Load your empirical crop harvest map raster (must be in EPSG:32630 metric coordinates)
harvest_map <-  rast("data/cassava_harvest_tonnes_2025.tif")

# Model Parameters
damage_coefficient <- 0.65  # Assumes a 65% aggregate yield loss penalty on infected hosts
price_per_tonne_usd <- 120.00 # Target farmgate price per metric ton of raw cassava tuberous roots

#-----------------------------------------------------------------------------
# 2. Extract Harvest Baseline to Simulation Nodes
#-----------------------------------------------------------------------------
# Convert animation_df back to an active spatial sf collection to sample pixels
nodes_for_extraction <- animation_df %>% 
  filter(Year == 1) %>% 
  st_as_sf(coords = c("x", "y"), crs = 32630) %>% 
  vect() # Convert to SpatVector for speed in terra

# Sample the exact underlying harvest metric for each specific farm node coordinate
extracted_harvest <- terra::extract(harvest_map, nodes_for_extraction)

# Join the baseline harvest capacities back to our main simulation tracking data frame
node_harvest_lookup <- tibble(
  node_id = 1:nrow(nodes_for_extraction),
  potential_harvest_tonnes = extracted_harvest[, 2] # Drops ID column, isolates raster values
) %>% 
  mutate(potential_harvest_tonnes = replace_na(potential_harvest_tonnes, 0))

#-----------------------------------------------------------------------------
# 3. Compute Annual Spatiotemporal Damage Metrics
#-----------------------------------------------------------------------------
cat(">>> Processing multi-year macroeconomic damage matrices...\n")

loss_analysis_df <- animation_df %>%
  mutate(node_id = rep(1:N_nodes, total_years)) %>% 
  left_join(node_harvest_lookup, by = "node_id") %>% 
  mutate(
    # Compute physical harvest reduction (tonnes)
    harvest_loss_tonnes = if_else(Status == "Infected", 
                                  potential_harvest_tonnes * damage_coefficient, 
                                  0),
    # Translate physical reduction into local market currency impacts
    economic_loss_usd = harvest_loss_tonnes * price_per_tonne_usd
  )

#-----------------------------------------------------------------------------
# 4. Generate Country-Wide Annual Ledger Reports
#-----------------------------------------------------------------------------
annual_summary_ledger <- loss_analysis_df %>% 
  group_by(Year) %>% 
  summarise(
    total_active_epicenters = sum(Status == "Infected"),
    cumulative_harvest_lost_tonnes = sum(harvest_loss_tonnes),
    cumulative_economic_loss_usd = sum(economic_loss_usd),
    .groups = "drop"
  )

# Save the dataset to outputs directory
write_csv(annual_summary_ledger, "outputs/epidemic_annual_economic_impact_ledger.csv")
print(head(annual_summary_ledger, 10))

#-----------------------------------------------------------------------------
# 5. Optional: Plot the Macroscopic Loss Curve
#-----------------------------------------------------------------------------
loss_curve_plot <- ggplot(annual_summary_ledger, aes(x = Year, y = cumulative_economic_loss_usd / 1e6)) +
  geom_line(color = "#cb181d", size = 1.2) +
  geom_point(color = "#1a1a1a", size = 2) +
  theme_minimal() +
  labs(
    title = "West Africa: Simulated Cumulative Epidemic Loss Trajectory",
    subtitle = paste0("Assumes a ", damage_coefficient*100, "% Yield Cut | $", price_per_tonne_usd, "/Tonne Root Pricing"),
    x = "Simulation Timeline (Years)",
    y = "Loss Scope (Millions of USD)"
  ) +
  theme(plot.title = element_text(face="bold", size=12))

ggsave("outputs/simulated_economic_loss_curve.png", plot = loss_curve_plot, width = 8, height = 5, dpi = 300)