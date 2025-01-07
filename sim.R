library(tidyverse)
library(mvtnorm)
library(norm)
library(mice)
library(naniar)
library(pbapply)

source("sim_functions.R")

### EVALUATION OF MEAN, COV, and COR

# iterated params
params_grid <- expand.grid(rho = seq(0,.3,.7),
                           type = c("MCAR", "MAR", "MNAR"))
# set params
B <- 1000 # sims per scenario
n <- 1000 # rows
d <- 3 # cols
mu <- 0 # colmeans
sigma <- 1 # sigma
size <- .10 # % missing, or percentile missing if MAR/MNAR

results <- data.frame(
  rho = numeric(),
  type = character(),
  res = I(list()),
  test = I(list())
)

for (row in 1:nrow(params_grid)) {
  rho <- params_grid$rho[row]
  type <- params_grid$type[row]
  MAR_rule <- get_mar_rule(type, size)
  
  # Generate parameters
  gen_params <- list(n = n, d = d, mu = mu, sigma = sigma, rho = rho)
  na_params <- list(size = size, type = type, MAR_rule = MAR_rule)
  
  # Run simulations
  sim <- pbreplicate(B, impute_sim(gen_params, na_params), simplify = FALSE)
  
  # Process results
  res <- process_results(sim, methods, statistics)
  test <- extract_mcar_values(sim)
  
  # Save to results data frame
  results <- results %>%
    add_row(rho = rho, type = type, res = list(res), test = list(test))
}

full_results <- combine_results(results)
full_results$means
full_results$cov
full_results$cor

# Calculate column statistics
cov_stats <- calculate_statistics(full_results$cov)
cor_stats <- calculate_statistics(full_results$cor)

# saveRDS(full_results, file = "simulation_results.RDS")







