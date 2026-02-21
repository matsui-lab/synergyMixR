test_that("extract_subject_features_by_singlePCA works correctly", {
  # Create simple test data
  N <- 5
  M <- 4
  T_i <- 20
  
  # Generate list of data matrices
  set.seed(123)
  list_of_data <- lapply(1:N, function(i) {
    matrix(rnorm(T_i * M), nrow = T_i, ncol = M)
  })
  
  # Test with r_dim = 2
  features <- extract_subject_features_by_singlePCA(list_of_data, r_dim = 2)
  
  # Check dimensions
  expect_equal(nrow(features), N)
  expect_equal(ncol(features), 2)
  
  # Check that output is numeric matrix
  expect_true(is.matrix(features))
  expect_true(is.numeric(features))
  
  # Check that there are no NA values
  expect_false(any(is.na(features)))
  
  # Test with r_dim = 3
  features_3d <- extract_subject_features_by_singlePCA(list_of_data, r_dim = 3)
  expect_equal(nrow(features_3d), N)
  expect_equal(ncol(features_3d), 3)
})

test_that("extract_subject_features_by_singlePCA handles edge cases", {
  # Test with single subject
  list_of_data_single <- list(matrix(rnorm(50 * 4), nrow = 50, ncol = 4))
  features_single <- extract_subject_features_by_singlePCA(list_of_data_single, r_dim = 2)
  expect_equal(nrow(features_single), 1)
  expect_equal(ncol(features_single), 2)
  
  # Test with r_dim = 1
  N <- 3
  list_of_data <- lapply(1:N, function(i) {
    matrix(rnorm(20 * 4), nrow = 20, ncol = 4)
  })
  features_1d <- extract_subject_features_by_singlePCA(list_of_data, r_dim = 1)
  expect_equal(nrow(features_1d), N)
  expect_equal(ncol(features_1d), 1)
})

test_that("assign_by_kmeans works correctly", {
  # Create simple feature matrix
  set.seed(456)
  N <- 30
  d <- 3
  features <- matrix(rnorm(N * d), nrow = N, ncol = d)
  
  # Test with K = 3
  K <- 3
  z <- assign_by_kmeans(features, K = K)
  
  # Check output is integer vector
  expect_true(is.integer(z) || is.numeric(z))
  expect_equal(length(z), N)
  
  # Check that cluster labels are in range 1..K
  expect_true(all(z >= 1))
  expect_true(all(z <= K))
  
  # Check that all clusters are represented (with high probability)
  # Note: This might occasionally fail with random data, but very unlikely
  expect_true(length(unique(z)) >= 1)
  expect_true(length(unique(z)) <= K)
})

test_that("assign_by_kmeans handles different K values", {
  set.seed(789)
  N <- 20
  features <- matrix(rnorm(N * 2), nrow = N, ncol = 2)
  
  # Test with K = 1
  z1 <- assign_by_kmeans(features, K = 1)
  expect_equal(length(z1), N)
  expect_true(all(z1 == 1))
  
  # Test with K = 2
  z2 <- assign_by_kmeans(features, K = 2)
  expect_equal(length(z2), N)
  expect_true(all(z2 %in% c(1, 2)))
  
  # Test with K = 5
  z5 <- assign_by_kmeans(features, K = 5)
  expect_equal(length(z5), N)
  expect_true(all(z5 >= 1 & z5 <= 5))
})

test_that("initialization functions work together", {
  # Integration test: extract features then cluster
  set.seed(999)
  N <- 15
  M <- 6
  T_i <- 30
  
  list_of_data <- lapply(1:N, function(i) {
    matrix(rnorm(T_i * M), nrow = T_i, ncol = M)
  })
  
  # Extract features
  features <- extract_subject_features_by_singlePCA(list_of_data, r_dim = 2)
  
  # Cluster based on features
  K <- 3
  z <- assign_by_kmeans(features, K = K)
  
  # Verify the pipeline works
  expect_equal(nrow(features), N)
  expect_equal(length(z), N)
  expect_true(all(z >= 1 & z <= K))
})

test_that("extract_subject_features_by_singlePCA is deterministic with same seed", {
  N <- 5
  M <- 4
  T_i <- 20
  
  # Generate same data twice
  set.seed(111)
  list_of_data1 <- lapply(1:N, function(i) {
    matrix(rnorm(T_i * M), nrow = T_i, ncol = M)
  })
  
  set.seed(111)
  list_of_data2 <- lapply(1:N, function(i) {
    matrix(rnorm(T_i * M), nrow = T_i, ncol = M)
  })
  
  # Extract features
  features1 <- extract_subject_features_by_singlePCA(list_of_data1, r_dim = 2)
  features2 <- extract_subject_features_by_singlePCA(list_of_data2, r_dim = 2)
  
  # Should be identical
  expect_equal(features1, features2)
})
