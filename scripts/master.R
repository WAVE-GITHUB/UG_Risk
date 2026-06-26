#==================================================
# Master file to run the entire piepline
#==================================================

source("scripts/data_cleaning.R")
source("scripts/powerLaw.R")
source("scripts/cassMap.R")
source("scripts/roads.R")
source("scripts/riskMap.R")
source("scripts/network_nodes.R")
source("scripts/ibm_smc_abc.R")
# NEW: Generate and render the 30-year spatiotemporal projection animation
source("scripts/simulate_progression.R")
source("scripts/compute_economic_losses.R")  