# Generate Test Data Fixtures for synergyMixR
# This script creates standardized test datasets for unit and integration testing

#' Generate Small Test Dataset
#' 
#' Creates a minimal dataset for unit testing with known properties
#' @param N Number of subjects (default: 5)
#' @param T_i Number of time points per subject (default: 20)
#' @param M Number of muscles/channels (default: 4)
#' @param K True number of clusters (default: 2)
#' @param r True number of factors (default: 2)
#' @param seed Random seed for reproducibility (default: 123)
#' @return List of data matrices and true parameters
generate_small_test_data <- function(N = 5, T_i = 20, M = 4, K = 2, r = 2, seed = 123) {
  set.seed(seed)
  
  # True cluster assignments
  z_true <- rep(1:K, length.out = N)
  
  # True parameters for each cluster
  mu_true <- list()
  Lambda_true <- list()
  Psi_true <- list()
  
  for (k in 1:K) {
    mu_true[[k]] <- rnorm(M, mean = k, sd = 0.5)
    Lambda_true[[k]] <- matrix(rnorm(M * r, mean = 0, sd = 1), M, r)
    Psi_true[[k]] <- diag(runif(M, 0.1, 0.5))
  }
  
  # Generate data
  list_of_data <- list()
  for (i in 1:N) {
    k <- z_true[i]
    
    # Generate factor scores
    Z_i <- matrix(rnorm(T_i * r), T_i, r)
    
    # Generate observations
    X_i <- Z_i %*% t(Lambda_true[[k]]) + 
           matrix(rep(mu_true[[k]], T_i), T_i, M, byrow = TRUE) +
           matrix(rnorm(T_i * M, 0, sqrt(diag(Psi_true[[k]]))), T_i, M)
    
    list_of_data[[i]] <- X_i
  }
  
  list(
    data = list_of_data,
    true_params = list(
      z = z_true,
      mu = mu_true,
      Lambda = Lambda_true,
      Psi = Psi_true,
      K = K,
      r = r,
      N = N,
      T_i = T_i,
      M = M
    )
  )
}

#' Generate Medium Test Dataset
#' 
#' Creates a medium-sized dataset for integration testing
#' @param N Number of subjects (default: 20)
#' @param T_i Number of time points per subject (default: 50)
#' @param M Number of muscles/channels (default: 8)
#' @param K True number of clusters (default: 3)
#' @param r True number of factors (default: 3)
#' @param seed Random seed for reproducibility (default: 456)
#' @return List of data matrices and true parameters
generate_medium_test_data <- function(N = 20, T_i = 50, M = 8, K = 3, r = 3, seed = 456) {
  set.seed(seed)
  
  # More realistic cluster assignments (unbalanced)
  cluster_probs <- c(0.4, 0.35, 0.25)
  z_true <- sample(1:K, N, replace = TRUE, prob = cluster_probs)
  
  # True parameters with more realistic structure
  mu_true <- list()
  Lambda_true <- list()
  Psi_true <- list()
  
  for (k in 1:K) {
    # Cluster-specific means
    mu_true[[k]] <- rnorm(M, mean = k * 0.5, sd = 0.3)
    
    # Structured loadings (some muscles more active in certain factors)
    Lambda_k <- matrix(0, M, r)
    for (j in 1:r) {
      # Each factor loads on a subset of muscles
      active_muscles <- sample(1:M, size = ceiling(M * 0.6))
      Lambda_k[active_muscles, j] <- rnorm(length(active_muscles), 
                                          mean = ifelse(k == j, 1.5, 0.5), 
                                          sd = 0.3)
    }
    Lambda_true[[k]] <- Lambda_k
    
    # Heteroscedastic noise
    Psi_true[[k]] <- diag(runif(M, 0.1, 0.8))
  }
  
  # Generate data with temporal correlation
  list_of_data <- list()
  for (i in 1:N) {
    k <- z_true[i]
    
    # Generate factor scores with temporal smoothness
    Z_i <- matrix(0, T_i, r)
    Z_i[1, ] <- rnorm(r)
    for (t in 2:T_i) {
      Z_i[t, ] <- 0.7 * Z_i[t-1, ] + 0.3 * rnorm(r)
    }
    
    # Generate observations
    X_i <- Z_i %*% t(Lambda_true[[k]]) + 
           matrix(rep(mu_true[[k]], T_i), T_i, M, byrow = TRUE) +
           matrix(rnorm(T_i * M, 0, sqrt(diag(Psi_true[[k]]))), T_i, M)
    
    list_of_data[[i]] <- X_i
  }
  
  list(
    data = list_of_data,
    true_params = list(
      z = z_true,
      mu = mu_true,
      Lambda = Lambda_true,
      Psi = Psi_true,
      K = K,
      r = r,
      N = N,
      T_i = T_i,
      M = M,
      cluster_probs = cluster_probs
    )
  )
}

#' Generate Large Test Dataset
#' 
#' Creates a large dataset for performance testing
#' @param N Number of subjects (default: 100)
#' @param T_i Number of time points per subject (default: 100)
#' @param M Number of muscles/channels (default: 12)
#' @param K True number of clusters (default: 4)
#' @param r True number of factors (default: 4)
#' @param seed Random seed for reproducibility (default: 789)
#' @return List of data matrices and true parameters
generate_large_test_data <- function(N = 100, T_i = 100, M = 12, K = 4, r = 4, seed = 789) {
  set.seed(seed)
  
  # Realistic cluster distribution
  cluster_probs <- c(0.3, 0.3, 0.25, 0.15)
  z_true <- sample(1:K, N, replace = TRUE, prob = cluster_probs)
  
  # Complex parameter structure
  mu_true <- list()
  Lambda_true <- list()
  Psi_true <- list()
  
  for (k in 1:K) {
    # Cluster-specific baseline activation
    mu_true[[k]] <- abs(rnorm(M, mean = k * 0.3, sd = 0.2))
    
    # Hierarchical factor structure
    Lambda_k <- matrix(0, M, r)
    
    # Global factor (affects all muscles)
    Lambda_k[, 1] <- rnorm(M, mean = 0.5, sd = 0.2)
    
    # Cluster-specific factors
    for (j in 2:r) {
      muscle_groups <- split(1:M, ceiling((1:M) / (M/r)))
      group_idx <- ((j-2) %% length(muscle_groups)) + 1
      active_muscles <- muscle_groups[[group_idx]]
      
      Lambda_k[active_muscles, j] <- rnorm(length(active_muscles), 
                                          mean = ifelse(k == j-1, 1.2, 0.3), 
                                          sd = 0.25)
    }
    Lambda_true[[k]] <- Lambda_k
    
    # Realistic noise structure
    Psi_true[[k]] <- diag(rgamma(M, shape = 2, rate = 4))
  }
  
  # Generate data with realistic properties
  list_of_data <- list()
  for (i in 1:N) {
    k <- z_true[i]
    
    # Generate factor scores with AR(1) structure
    Z_i <- matrix(0, T_i, r)
    Z_i[1, ] <- rnorm(r, sd = 0.5)
    for (t in 2:T_i) {
      Z_i[t, ] <- 0.8 * Z_i[t-1, ] + sqrt(1 - 0.8^2) * rnorm(r, sd = 0.5)
    }
    
    # Add some non-stationarity
    trend <- seq(0, 0.1, length.out = T_i)
    Z_i <- Z_i + outer(trend, rep(1, r))
    
    # Generate observations
    X_i <- Z_i %*% t(Lambda_true[[k]]) + 
           matrix(rep(mu_true[[k]], T_i), T_i, M, byrow = TRUE) +
           matrix(rnorm(T_i * M, 0, sqrt(diag(Psi_true[[k]]))), T_i, M)
    
    # Ensure non-negative values (realistic for EMG-like data)
    X_i <- pmax(X_i, 0.01)
    
    list_of_data[[i]] <- X_i
  }
  
  list(
    data = list_of_data,
    true_params = list(
      z = z_true,
      mu = mu_true,
      Lambda = Lambda_true,
      Psi = Psi_true,
      K = K,
      r = r,
      N = N,
      T_i = T_i,
      M = M,
      cluster_probs = cluster_probs
    )
  )
}

#' Generate Edge Case Test Data
#' 
#' Creates datasets for testing edge cases and boundary conditions
#' @return List of various edge case datasets
generate_edge_case_data <- function() {
  edge_cases <- list()
  
  # Single subject
  edge_cases$single_subject <- generate_small_test_data(N = 1, K = 1, seed = 111)
  
  # Single cluster
  edge_cases$single_cluster <- generate_small_test_data(K = 1, seed = 222)
  
  # Single factor
  edge_cases$single_factor <- generate_small_test_data(r = 1, seed = 333)
  
  # Very small time series
  edge_cases$short_timeseries <- generate_small_test_data(T_i = 5, seed = 444)
  
  # High dimensional
  edge_cases$high_dimensional <- generate_small_test_data(M = 20, r = 5, seed = 555)
  
  # Unbalanced clusters
  set.seed(666)
  N <- 10
  z_unbalanced <- c(rep(1, 8), rep(2, 2))  # Very unbalanced
  edge_cases$unbalanced <- generate_small_test_data(N = N, seed = 666)
  edge_cases$unbalanced$true_params$z <- z_unbalanced
  
  edge_cases
}

# Generate and save all test datasets
if (!interactive()) {
  # Create all test datasets
  small_data <- generate_small_test_data()
  medium_data <- generate_medium_test_data()
  large_data <- generate_large_test_data()
  edge_cases <- generate_edge_case_data()
  
  # Save to RDS files
  saveRDS(small_data, "fixtures/small_test_data.rds")
  saveRDS(medium_data, "fixtures/medium_test_data.rds")
  saveRDS(large_data, "fixtures/large_test_data.rds")
  saveRDS(edge_cases, "fixtures/edge_case_data.rds")
  
  cat("Test data fixtures generated successfully!\n")
  cat("Files created:\n")
  cat("- fixtures/small_test_data.rds\n")
  cat("- fixtures/medium_test_data.rds\n")
  cat("- fixtures/large_test_data.rds\n")
  cat("- fixtures/edge_case_data.rds\n")
}
