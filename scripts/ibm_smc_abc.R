#=============================================================================
# scripts/ibm_smc_abc.R
# Empirically-Grounded Network SI Model & SMC-ABC Calibration Engine
# Production-Grade Covariance-Insulated Edition
#=============================================================================

library(tidyverse)
library(EasyABC)

#-----------------------------------------------------------------------------
# 1. Verification and Workspace Initialization ----
#-----------------------------------------------------------------------------
cat("-> Verifying and initializing baseline workspace elements...\n")
if(!file.exists("outputs/nodes.rds") || !file.exists("outputs/movement_matrix.rds")) {
  stop("CRITICAL PIPELINE FAULT: Network nodes or border-weighted matrices are missing.")
}

raw_nodes_df        <- readRDS("outputs/nodes.rds")
raw_movement_matrix <- readRDS("outputs/movement_matrix.rds")

#-----------------------------------------------------------------------------
# 2. Define the Fully Insulated Network Simulation Engine ----
#-----------------------------------------------------------------------------
create_insulated_simulator <- function(nodes, move_mat) {
  
  # Clean, un-named matrix and setup variables for cross-environment safety
  env_matrix <- unname(as.matrix(move_mat))
  env_N      <- nrow(nodes)
  env_max_t  <- 30
  env_seeds  <- order(nodes$risk, decreasing = TRUE)[1:3]
  
  function(params) {
    # Extract the clean parameter scalar
    transmission_scale <- as.numeric(params[1])
    
    # 1L = Susceptible, 2L = Infected
    status <- rep(1L, env_N) 
    status[env_seeds] <- 2L
    
    infectious_curve    <- numeric(env_max_t)
    infectious_curve[1] <- length(env_seeds)
    
    # Time iteration execution loop
    for (t in 2:env_max_t) {
      current_infected    <- which(status == 2L)
      current_susceptible <- which(status == 1L)
      
      if (length(current_infected) == 0 || length(current_susceptible) == 0) {
        infectious_curve[seq(t, env_max_t)] <- length(current_infected)
        break
      }
      
      # Vectorized transmission matrix calculation
      sub_matrix <- env_matrix[current_susceptible, current_infected, drop = FALSE]
      
      # Joint probability of escaping infection across all active links
      escape_probabilities <- exp(rowSums(log(1 - pmin(sub_matrix * transmission_scale, 0.9999))))
      infection_risk       <- 1 - escape_probabilities
      
      # Stochastic state transition step
      new_infections <- current_susceptible[runif(length(current_susceptible)) < infection_risk]
      if (length(new_infections) > 0) {
        status[new_infections] <- 2L
      }
      
      infectious_curve[t] <- sum(status == 2L)
    }
    
    # CRITICAL COVARIANCE GUARDRAIL: 
    # Add a microscopic, parameter-dependent jitter to the tail of the output vector.
    # This breaks mathematical ties and guarantees non-zero variance across particles,
    # preventing EasyABC's internal covariance metrics from collapsing into NULL.
    jitter_factor <- (transmission_scale * 1e-8)
    clean_vector  <- as.numeric(infectious_curve) + jitter_factor
    
    return(clean_vector)
  }
}

# Instantiate the environment-safe simulator
simulate_empirical_si_network <- create_insulated_simulator(raw_nodes_df, raw_movement_matrix)

#-----------------------------------------------------------------------------
# 3. Handle Observed Target Validation Data ----
#-----------------------------------------------------------------------------
cat("-> Generating validation target trajectories...\n")
true_calibration_scale <- c(0.45) 
observed_trajectory    <- simulate_empirical_si_network(true_calibration_scale)

#-----------------------------------------------------------------------------
# 4. Execute the Lenormand Adaptive SMC-ABC Optimization Sequence ----
#-----------------------------------------------------------------------------
prior_distribution <- list(
  c("unif", 0.001, 5.0)
)

cat(">>> Launching Sequential Monte Carlo ABC across insulated environment grids...\n")

abc_smc_results <- ABC_sequential(
  method              = "Lenormand",
  model               = simulate_empirical_si_network,
  prior               = prior_distribution,
  summary_stat_target = observed_trajectory,
  nb_simul            = 250,            # Particle evaluation population size
  alpha               = 0.2,            # Top 20% tolerance selector
  verbose             = FALSE
)

#-----------------------------------------------------------------------------
# 5. Extract and Save Calibrated Posteriors ----
#-----------------------------------------------------------------------------
cat("-> Exporting calibrated network parameters...\n")

posterior_weights      <- abc_smc_results$weights
posterior_param_matrix <- as.matrix(abc_smc_results$param)

optimal_scaler <- sum(posterior_param_matrix[, 1] * posterior_weights) / sum(posterior_weights)

optimal_parameters <- tibble(
  parameter        = "transmission_scale",
  optimal_value    = optimal_scaler,
  calibration_date = Sys.Date()
)

saveRDS(optimal_parameters, "outputs/optimal_abc_parameters.rds")
write_csv(optimal_parameters, "outputs/optimal_abc_parameters.csv")

cat("✔ Pipeline calibration finalized successfully.\n\n")
print(optimal_parameters)