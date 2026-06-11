#=============================================================================
# scripts/ibm_smc_abc.R
# Empirically-Grounded Individual-Based SI Model & SMC-ABC Calibration
#=============================================================================

library(tidyverse)
library(EasyABC)

#-----------------------------------------------------------------------------
# 1. Load Real Spatial Data Derived from Master Pipeline ----
#-----------------------------------------------------------------------------
if(!file.exists("outputs/nodes.rds") || !file.exists("outputs/distance_matrix.rds")) {
  stop("Missing structural inputs. Please run 'scripts/network_nodes.R' first.")
}

empirical_nodes <- readRDS("outputs/nodes.rds")
empirical_dist_matrix <- readRDS("outputs/distance_matrix.rds")

# Dynamically set simulation dimensions to match your landscape network layers
N_nodes <- nrow(empirical_nodes)
max_time <- 30  # Longitudinal tracking window iterations

#-----------------------------------------------------------------------------
# 2. Define the Empirically-Grounded Individual-Based Model ----
#-----------------------------------------------------------------------------
simulate_empirical_si_ibm <- function(params) {
  beta <- params[1]
  gamma_dist <- params[2]
  
  # Initialize node status based on real landscape shape
  # Setting nodes with the highest risk baseline as initial index infection sites
  status <- rep("S", N_nodes)
  initial_infected_indices <- order(empirical_nodes$risk, decreasing = TRUE)[1:3]
  status[initial_infected_indices] <- "I"
  
  # Vector to track cumulative epidemic curve trajectory
  infectious_curve <- numeric(max_time)
  infectious_curve[1] <- length(initial_infected_indices)
  
  # Temporal simulation loop
  for (t in 2:max_time) {
    current_infected <- which(status == "I")
    current_susceptible <- which(status == "S")
    
    if (length(current_infected) == 0 || length(current_susceptible) == 0) {
      infectious_curve[t:max_time] <- length(current_infected)
      break
    }
    
    # Force of infection using the real distance matrix
    escape_probabilities <- sapply(current_susceptible, function(s_idx) {
      # Pull exact metric distances from distance matrix output
      distances <- empirical_dist_matrix[s_idx, current_infected]
      
      # Spatial Power-Law Kernel mediated by distance parameters
      kernel <- beta / (1 + (distances / 1000)^gamma_dist) # Distances scaled to km for stability
      prod(1 - kernel)
    })
    
    infection_risk <- 1 - escape_probabilities
    
    # Stochastic infection state transitions
    new_infections <- current_susceptible[runif(length(current_susceptible)) < infection_risk]
    
    if (length(new_infections) > 0) {
      status[new_infections] <- "I"
    }
    
    infectious_curve[t] <- sum(status == "I")
  }
  
  return(infectious_curve)
}

#-----------------------------------------------------------------------------
# 3. Handle Observed Target Validation Data ----
#-----------------------------------------------------------------------------
# NOTE: Replace 'synthetic_target' with your real-world field outbreak timelines
# (e.g., aggregate observations from dataEpi over time)
# For runtime validity, we generate a target using representative baseline metrics:
true_baseline_params <- c(0.35, 1.6)
observed_trajectory <- simulate_empirical_si_ibm(true_baseline_params)

#-----------------------------------------------------------------------------
# 4. Execute the SMC-ABC Optimization Sequence (Lenormand Adaptive)
#-----------------------------------------------------------------------------
prior_distribution <- list(
  c("unif", 0.01, 2.0),  # beta bounds
  c("unif", 0.5, 4.0)    # gamma bounds
)

cat(">>> Launching Sequential Monte Carlo ABC across empirical grids...\n")

abc_smc_results <- ABC_sequential(
  method = "Lenormand",                  # Adaptive framework
  model = simulate_empirical_si_ibm,
  prior = prior_distribution,
  summary_stat_target = observed_trajectory,
  nb_simul = 100,                        
  alpha = 0.2,                           # Keeps the top 20% of particles
  verbose = FALSE
)

#-----------------------------------------------------------------------------
# 5. Extract and Save Calibrated Posteriors ----
#-----------------------------------------------------------------------------
posterior_weights <- abc_smc_results$weights
posterior_param_matrix <- abc_smc_results$param

# Compute weighted mean point estimates for optimal fit
optimal_beta  <- sum(posterior_param_matrix[, 1] * posterior_weights) / sum(posterior_weights)
optimal_gamma <- sum(posterior_param_matrix[, 2] * posterior_weights) / sum(posterior_weights)

optimal_parameters <- tibble(
  parameter = c("beta", "gamma_dist"),
  optimal_value = c(optimal_beta, optimal_gamma),
  calibration_date = Sys.Date()
)

# Export parameters to outputs directory for network execution scripts
saveRDS(optimal_parameters, "outputs/optimal_abc_parameters.rds")
write_csv(optimal_parameters, "outputs/optimal_abc_parameters.csv")

cat("✔ Calibration finalized. Optimal parameters saved to 'outputs/optimal_abc_parameters.rds'.\n\n")
print(optimal_parameters)