#!/usr/bin/env Rscript
# =============================================================================
# Monte Carlo Simulation: Welch's Two-Sample Test on the Counterexample Family
#
# Investigates the finite-sample rejection rate of Welch's t-test under the
# symmetric triangular-array / fixed family construction:
#   Group X: rare amplitude sqrt(m/8) with prob 8/(2m)
#   Group Y: rare amplitude sqrt(m/8) with prob 12/(2m)
#
# Reproducibility Seed: 110
# =============================================================================

set.seed(110)

# Function to run Welch simulation for a specific sample size m
simulate_welch <- function(m, n_sim = 200000, alpha = 0.05) {
  lambda_X <- 8
  lambda_Y <- 12
  
  p_rare_X <- lambda_X / (2 * m)
  p_rare_Y <- lambda_Y / (2 * m)
  
  rare_val <- sqrt(m / 8)
  bg_val   <- 1 / m
  
  rejections <- integer(n_sim)
  zero_var_count <- 0
  
  cat(sprintf("\nRunning simulation for sample size m = %d (Replications = %d)...\n", m, n_sim))
  
  # Chunked simulation to optimize memory and speed
  chunk_size <- 20000
  n_chunks <- ceiling(n_sim / chunk_size)
  
  total_done <- 0
  
  for (chunk in 1:n_chunks) {
    cur_n <- min(chunk_size, n_sim - total_done)
    if (cur_n <= 0) break
    
    # Generate Group X
    # Discrete choices: +rare, -rare, +bg, -bg
    # Probabilities: [p_rare_X, p_rare_X, (1-2*p_rare_X)/2, (1-2*p_rare_X)/2]
    prob_X <- c(p_rare_X, p_rare_X, (1 - 2*p_rare_X)/2, (1 - 2*p_rare_X)/2)
    vals_X <- c(rare_val, -rare_val, bg_val, -bg_val)
    
    X_mat <- matrix(sample(vals_X, size = cur_n * m, replace = TRUE, prob = prob_X),
                    nrow = cur_n, ncol = m)
    
    # Generate Group Y
    prob_Y <- c(p_rare_Y, p_rare_Y, (1 - 2*p_rare_Y)/2, (1 - 2*p_rare_Y)/2)
    vals_Y <- c(rare_val, -rare_val, bg_val, -bg_val)
    
    Y_mat <- matrix(sample(vals_Y, size = cur_n * m, replace = TRUE, prob = prob_Y),
                    nrow = cur_n, ncol = m)
    
    # Compute row-wise means
    mean_X <- rowMeans(X_mat)
    mean_Y <- rowMeans(Y_mat)
    
    # Compute row-wise sample variances (unbiased)
    # Fast row-variance calculation
    var_X <- rowSums((X_mat - mean_X)^2) / (m - 1)
    var_Y <- rowSums((Y_mat - mean_Y)^2) / (m - 1)
    
    denom_sq <- var_X / m + var_Y / m
    regular_domain <- (denom_sq > 0)
    
    # Welch degrees of freedom
    # nu = (var_X/m + var_Y/m)^2 / [ (var_X/m)^2/(m-1) + (var_Y/m)^2/(m-1) ]
    nu <- (denom_sq^2) / ( (var_X / m)^2 / (m - 1) + (var_Y / m)^2 / (m - 1) )
    
    # Welch t-statistic
    t_stat <- rep(0, cur_n)
    t_stat[regular_domain] <- (mean_X[regular_domain] - mean_Y[regular_domain]) / sqrt(denom_sq[regular_domain])
    
    # Two-sided critical values
    crit_val <- qt(1 - alpha / 2, df = nu)
    
    # Rejection decision rule (strictly 0 on degenerate domain)
    rejected <- regular_domain & (abs(t_stat) > crit_val)
    
    rejections[(total_done + 1):(total_done + cur_n)] <- as.integer(rejected)
    zero_var_count <- zero_var_count + sum(!regular_domain)
    
    total_done <- total_done + cur_n
  }
  
  emp_rate <- mean(rejections)
  se_rate  <- sqrt(emp_rate * (1 - emp_rate) / n_sim)
  ci_lower <- emp_rate - 1.96 * se_rate
  ci_upper <- emp_rate + 1.96 * se_rate
  
  cat(sprintf("  Results for m = %d:\n", m))
  cat(sprintf("    Empirical Rejection Rate: %.5f\n", emp_rate))
  cat(sprintf("    95%% Confidence Interval:   [%.5f, %.5f]\n", ci_lower, ci_upper))
  cat(sprintf("    Zero-variance events:     %d\n", zero_var_count))
  cat(sprintf("    Theoretical Lower Bound:  0.05191\n"))
  
  return(data.frame(m = m, emp_rate = emp_rate, se = se_rate, ci_lo = ci_lower, ci_hi = ci_upper))
}

# =============================================================================
# Execution Grid
# =============================================================================
cat("======================================================================\n")
cat("Starting Welch Test Finite-Truncation Simulation Study (seed = 110)\n")
cat("======================================================================\n")

sample_sizes <- c(50, 200, 1000, 5000)
results_list <- list()

for (i in seq_along(sample_sizes)) {
  results_list[[i]] <- simulate_welch(m = sample_sizes[i], n_sim = 200000)
}

final_results <- do.call(rbind, results_list)

cat("\n======================================================================\n")
cat("SUMMARY TABLE:\n")
cat("======================================================================\n")
print(final_results, row.names = FALSE)
cat("======================================================================\n")
cat("The empirical rejection rates are compared with the certified asymptotic lower bound 0.051913.\n")
cat("The lower bound is not asserted to equal the exact limiting rejection probability.\n")