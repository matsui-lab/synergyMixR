# Testing Guide for synergyMixR

## Overview

This document provides comprehensive guidance on testing the synergyMixR R package. The testing infrastructure is designed to ensure code quality, mathematical correctness, and reproducibility across different environments.

## Quick Start

### Running All Tests

```r
# From R console
devtools::test()

# Or from command line
R CMD check .
```

### Running Specific Test Files

```r
# Test specific functionality
testthat::test_file("tests/testthat/test-mfa-em-fit.R")
testthat::test_file("tests/testthat/test-mpca-em-fit.R")
testthat::test_file("tests/testthat/test-utility.R")
```

### Running Tests with Coverage

```r
# Generate coverage report
covr::package_coverage()

# View interactive coverage report
covr::report()
```

## Test Organization

### Test Structure

```
tests/
├── testthat.R                              # Test runner (required by R CMD check)
├── TESTING_STRATEGY.md                     # Comprehensive testing strategy
├── TESTING_GUIDE.md                        # This file
└── testthat/
    ├── fixtures/                           # Test data
    │   ├── generate_test_data.R           # Data generation script
    │   ├── small_test_data.rds            # Small dataset (N=5, T=20, M=4)
    │   ├── medium_test_data.rds           # Medium dataset (N=20, T=50, M=8)
    │   ├── large_test_data.rds            # Large dataset (N=100, T=100, M=12)
    │   └── edge_case_data.rds             # Edge cases (single cluster, etc.)
    ├── test-initialization_utils.R        # Initialization utilities tests
    ├── test-mfa-em-fit.R                  # MFA algorithm tests
    ├── test-mpca-em-fit.R                 # MPCA algorithm tests
    ├── test-utility.R                     # Utility function tests
    ├── test-utility-extended.R            # Extended utility tests (edge cases, numerical precision)
    ├── test-rotation.R                    # Factor rotation tests
    ├── test-flip.R                        # Sign alignment tests
    ├── test-best-fit.R                    # Model selection tests
    └── test-integration-workflows.R       # End-to-end workflow tests
```

### Test Categories

1. **Unit Tests**: Test individual functions in isolation
   - `test-mfa-em-fit.R`: Core MFA fitting functions
   - `test-mpca-em-fit.R`: Core MPCA fitting functions
   - `test-utility.R`: Utility and helper functions
   - `test-utility-extended.R`: Extended utility tests with edge cases and numerical precision
   - `test-rotation.R`: Factor rotation functions (varimax, score consistency)
   - `test-flip.R`: Sign alignment functions (Lambda and factor_scores correspondence)
   - `test-best-fit.R`: Model selection functions

2. **Integration Tests**: Test complete workflows
   - `test-integration-workflows.R`: End-to-end analysis pipelines

3. **Existing Tests**: Previously implemented tests
   - `test-initialization_utils.R`: Initialization and k-means utilities

## Test Data Fixtures

### Available Datasets

All test datasets are generated with known true parameters for validation:

1. **Small Test Data** (`small_test_data.rds`)
   - N=5 subjects, T=20 timepoints, M=4 muscles
   - K=2 clusters, r=2 factors
   - Use for: Quick unit tests, CI/CD

2. **Medium Test Data** (`medium_test_data.rds`)
   - N=20 subjects, T=50 timepoints, M=8 muscles
   - K=3 clusters, r=3 factors
   - Use for: Integration tests, realistic scenarios

3. **Large Test Data** (`large_test_data.rds`)
   - N=100 subjects, T=100 timepoints, M=12 muscles
   - K=4 clusters, r=4 factors
   - Use for: Performance tests, scalability validation

4. **Edge Case Data** (`edge_case_data.rds`)
   - Single subject, single cluster, single factor, etc.
   - Use for: Boundary condition testing

### Regenerating Test Data

```r
# Run the data generation script
source("tests/testthat/fixtures/generate_test_data.R")
```

## Writing New Tests

### Test Template

```r
# tests/testthat/test-new-feature.R

# Load test fixtures
small_data <- readRDS("fixtures/small_test_data.rds")

test_that("new_function works with valid inputs", {
  # Arrange
  list_of_data <- small_data$data
  true_params <- small_data$true_params
  
  # Act
  result <- new_function(list_of_data, param1 = value1)
  
  # Assert
  expect_type(result, "list")
  expect_true(is.finite(result$metric))
  expect_equal(length(result$output), true_params$N)
})

test_that("new_function handles edge cases", {
  edge_cases <- readRDS("fixtures/edge_case_data.rds")
  
  # Test with single cluster
  result <- new_function(edge_cases$single_cluster$data)
  expect_true(is.finite(result$metric))
})

test_that("new_function validates inputs", {
  # Test error handling
  expect_error(
    new_function(list()),
    "data"
  )
})
```

### Best Practices

1. **Use Descriptive Test Names**: Test names should clearly describe what is being tested
2. **Test One Thing**: Each test should focus on a single aspect of functionality
3. **Use Fixtures**: Load test data from fixtures rather than generating inline
4. **Test Edge Cases**: Include tests for boundary conditions and error cases
5. **Check Dimensions**: Verify output dimensions match expected values
6. **Validate Ranges**: Ensure outputs are within valid ranges (e.g., probabilities in [0,1])
7. **Test Determinism**: Verify reproducibility with fixed seeds

## Common Test Patterns

### Testing Model Fitting

```r
test_that("model_fit produces valid output structure", {
  fit <- model_fit(data, K = 2, r = 2)
  
  # Check structure
  expect_type(fit, "list")
  expect_named(fit, c("z", "pi", "mu", "Lambda", "Psi", "logLik"))
  
  # Check dimensions
  expect_equal(length(fit$z), N)
  expect_equal(length(fit$pi), K)
  
  # Check validity
  expect_true(all(fit$z %in% 1:K))
  expect_equal(sum(fit$pi), 1, tolerance = 1e-6)
  expect_true(is.finite(fit$logLik))
})
```

### Testing Mathematical Properties

```r
test_that("covariance matrix is positive definite", {
  Sigma <- compute_covariance(Lambda, Psi)
  
  # Check symmetry
  expect_equal(Sigma, t(Sigma), tolerance = 1e-10)
  
  # Check positive definiteness
  eigenvalues <- eigen(Sigma, only.values = TRUE)$values
  expect_true(all(eigenvalues > 0))
})
```

### Testing Convergence

```r
test_that("EM algorithm converges", {
  fit_short <- model_fit(data, max_iter = 5)
  fit_long <- model_fit(data, max_iter = 50)
  
  # Longer run should have equal or better log-likelihood
  expect_true(fit_long$logLik >= fit_short$logLik - 1e-6)
})
```

### Testing Reproducibility

```r
test_that("results are reproducible with fixed seed", {
  fit1 <- model_fit(data, seed = 123)
  fit2 <- model_fit(data, seed = 123)
  
  expect_equal(fit1$z, fit2$z)
  expect_equal(fit1$logLik, fit2$logLik, tolerance = 1e-10)
})
```

## Continuous Integration

### GitHub Actions Workflow

The package uses GitHub Actions for automated testing on multiple platforms:

- **Platforms**: Ubuntu (latest, devel, oldrel-1), macOS (latest), Windows (latest)
- **R Versions**: devel, release, oldrel-1
- **Checks**: R CMD check, lintr, test coverage

### Local Pre-commit Checks

Before committing, run:

```r
# Check package
devtools::check()

# Run tests
devtools::test()

# Check code style
lintr::lint_package()
```

## Debugging Failed Tests

### Running Tests Interactively

```r
# Load package in development mode
devtools::load_all()

# Load test fixtures
small_data <- readRDS("tests/testthat/fixtures/small_test_data.rds")

# Run code interactively
result <- mfa_em_fit(small_data$data, K = 2, r = 2, max_iter = 5)

# Inspect results
str(result)
summary(result)
```

### Common Issues

1. **Numerical Precision**: Use `tolerance` parameter in `expect_equal()`
2. **Random Initialization**: Set seeds for reproducibility
3. **Platform Differences**: Test on multiple platforms if possible
4. **Memory Issues**: Use smaller datasets for unit tests

## Performance Testing

### Benchmarking

```r
# Benchmark different implementations
microbenchmark::microbenchmark(
  mfa = mfa_em_fit(data, K = 2, r = 2, max_iter = 10),
  mpca = mixture_pca_em_fit(data, K = 2, r = 2, max_iter = 10),
  times = 10
)
```

### Profiling

```r
# Profile code execution
profvis::profvis({
  result <- mfa_em_fit(large_data$data, K = 3, r = 3, max_iter = 20)
})
```

## Test Coverage Goals

### Target Coverage

- **Overall Package**: ≥75% line coverage (target from Phase 2)
- **Core Algorithms**: ≥95% coverage (MFA/MPCA fitting)
- **Utilities**: ≥85% coverage
- **Post-processing Functions**: ≥75% coverage (rotation, flip, utility)
  - `rotation.R`: ≥75% (up from ~45%)
  - `flip.R`: ≥75% (up from ~30%)
  - `utility.R`: ≥75% (up from ~55%)
- **Visualization**: ≥70% coverage (harder to test)

### File-Specific Coverage Targets (Phase 2 - Issue #27)

The following files have enhanced test coverage as of Phase 2:

| File | Previous Coverage | Target Coverage | Test File |
|------|------------------|-----------------|-----------|
| `rotation.R` | ~45% | ≥75% | `test-rotation.R` |
| `flip.R` | ~30% | ≥75% | `test-flip.R` |
| `utility.R` | ~55% | ≥75% | `test-utility.R`, `test-utility-extended.R` |

### Coverage Policy

1. **New Code**: All new functions must have ≥90% test coverage
2. **Bug Fixes**: Add regression tests for all bug fixes
3. **Refactoring**: Maintain or improve coverage during refactoring
4. **Edge Cases**: Explicitly test boundary conditions and extreme values
5. **Numerical Precision**: Test functions with near-zero variance, large values, and small noise

### Checking Coverage

```r
# Generate coverage report
cov <- covr::package_coverage()

# View summary
cov

# View detailed report
covr::report(cov)

# Check specific files
covr::file_coverage("R/mfa_em_fit.R", "tests/testthat/test-mfa-em-fit.R")
covr::file_coverage("R/rotation.R", "tests/testthat/test-rotation.R")
covr::file_coverage("R/flip.R", "tests/testthat/test-flip.R")
covr::file_coverage("R/utility.R", c("tests/testthat/test-utility.R", 
                                      "tests/testthat/test-utility-extended.R"))
```

## Troubleshooting

### Tests Pass Locally but Fail in CI

1. Check platform-specific issues (Windows vs. Unix)
2. Verify all dependencies are installed
3. Check for hardcoded paths
4. Ensure tests don't depend on external resources

### Tests Are Too Slow

1. Use smaller datasets for unit tests
2. Reduce number of iterations in tests
3. Use `skip_on_cran()` for slow tests
4. Consider parallel test execution

### Flaky Tests

1. Set random seeds explicitly
2. Increase tolerance for numerical comparisons
3. Avoid timing-dependent assertions
4. Check for race conditions in parallel code

## Additional Resources

- [testthat documentation](https://testthat.r-lib.org/)
- [R Packages book - Testing chapter](https://r-pkgs.org/testing-basics.html)
- [covr package documentation](https://covr.r-lib.org/)
- [TESTING_STRATEGY.md](TESTING_STRATEGY.md) - Comprehensive testing strategy

## Getting Help

If you encounter issues with testing:

1. Check this guide and TESTING_STRATEGY.md
2. Review existing test files for examples
3. Run tests interactively to debug
4. Check GitHub Actions logs for CI failures
5. Open an issue on GitHub with reproducible example
