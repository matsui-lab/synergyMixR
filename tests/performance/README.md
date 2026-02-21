# Performance Tests for OpenMP Parallelization

## Overview

This directory contains performance tests to validate the OpenMP parallelization implementation in the `synergyMixR` package. These tests verify that:

1. **Reproducibility**: Results are identical across different thread counts when using the same seed
2. **Performance Scaling**: Execution time decreases (or stays reasonable) with increased thread count
3. **Data Size Scaling**: Performance benefits scale appropriately with data size
4. **Parameter Validation**: Thread count parameter is handled correctly
5. **Seed Consistency**: The seed parameter ensures reproducible results

## Running the Tests

### Quick Start

From the repository root:

```bash
Rscript tests/performance/test_parallel_performance.R
```

Or from within R:

```r
source("tests/performance/test_parallel_performance.R")
```

### Expected Output

The test suite will run 10 tests (5 test categories × 2 methods) and display:

```
========================================
Performance Test Suite for OpenMP
========================================

System Information:
  R version: R version 4.5.1 (2025-06-13)
  Platform: x86_64-pc-linux-gnu
  CPU cores: 8
  synergyMixR version: 0.3.0

Test 1: Reproducibility with Different Thread Counts
------------------------------------------------------
Testing MFA (N=50, K=2, r=2)...
[✓ PASS] MFA reproducibility
        logLik: 1-thread=-1234.567, 2-threads=-1234.567, 4-threads=-1234.567 (max diff=1.23e-12)
Testing MPCA (N=50, K=2, r=2)...
[✓ PASS] MPCA reproducibility
        logLik: 1-thread=-1234.567, 2-threads=-1234.567, 4-threads=-1234.567 (max diff=1.23e-12)

Test 2: Performance Scaling with Thread Count
------------------------------------------------------
Testing MFA performance scaling (N=100, K=2, r=2)...
  1 thread : 2.47 s
  2 threads: 2.15 s
  4 threads: 2.21 s
  Speedup: 2-threads=1.15x, 4-threads=1.12x
[✓ PASS] MFA performance scaling
        Speedup: 2-threads=1.15x, 4-threads=1.12x (baseline: 2.47 s)
...

========================================
Test Summary
========================================
Total tests: 10
Passed: 10
Failed: 0

✓ All performance tests passed!
```

### Exit Codes

- **0**: All tests passed
- **1**: One or more tests failed

## Test Descriptions

### Test 1: Reproducibility with Different Thread Counts

**Purpose**: Verify that using the same `cpp_seed` produces identical results regardless of thread count.

**Method**: 
- Generate synthetic data with fixed seed
- Run `mfa_em_fit()` and `mixture_pca_em_fit()` with 1, 2, and 4 threads
- Compare log-likelihood values (should be identical within numerical tolerance)

**Pass Criteria**: Maximum difference < 1e-10

### Test 2: Performance Scaling with Thread Count

**Purpose**: Verify that increasing thread count provides speedup (or at least doesn't significantly degrade performance).

**Method**:
- Generate synthetic data (N=100, K=2, r=2)
- Measure execution time with 1, 2, and 4 threads
- Calculate speedup relative to single-threaded execution

**Pass Criteria**: Speedup ≥ 0.8× (allows for modest speedup due to overhead)

**Note**: Modest speedup (0.9-1.15×) is expected due to:
- Sequential M-step
- Memory bandwidth limitations
- Thread synchronization overhead
- Amdahl's Law

### Test 3: Performance Scaling with Data Size

**Purpose**: Verify that parallelization benefits scale appropriately with data size.

**Method**:
- Test with N = 20, 50, 100 subjects
- Measure execution time with 1 and 4 threads for each data size
- Calculate average speedup across data sizes

**Pass Criteria**: Average speedup ≥ 0.8×

**Expected Behavior**: Larger datasets should show equal or better speedup compared to smaller datasets.

### Test 4: Thread Count Parameter Validation

**Purpose**: Verify that the `n_threads` parameter is handled correctly for various values.

**Method**:
- Test with n_threads = 1, 2, 4, 8
- Verify that each configuration produces valid results
- Check that log-likelihood is finite and non-null

**Pass Criteria**: All thread counts produce valid results without errors

### Test 5: Seed Parameter Consistency

**Purpose**: Verify that using the same seed produces identical results across multiple runs.

**Method**:
- Run the same fit 3 times with identical parameters and seed
- Compare log-likelihood values

**Pass Criteria**: All runs produce identical results (max difference < 1e-10)

## Interpreting Results

### All Tests Pass

If all tests pass, the OpenMP parallelization is working correctly:
- ✓ Results are reproducible
- ✓ Performance scales reasonably with thread count
- ✓ Thread parameter handling is robust
- ✓ Seed parameter ensures consistency

### Some Tests Fail

If tests fail, investigate:

1. **Reproducibility failures**: 
   - Check OpenMP installation
   - Verify RNG seeding implementation
   - Review thread-safety of shared data structures

2. **Performance scaling failures**:
   - Check system has multiple CPU cores available
   - Verify OpenMP is enabled (`OMP_NUM_THREADS` environment variable)
   - Review system load (other processes may affect performance)
   - Note: Speedup < 0.8× may indicate overhead dominates for small workloads

3. **Parameter validation failures**:
   - Check for segmentation faults or crashes
   - Review error messages for clues
   - Verify OpenMP configuration in `src/Makevars`

4. **Seed consistency failures**:
   - Check RNG implementation in C++ code
   - Verify seed is properly passed to C++ functions
   - Review thread-local RNG initialization

## Customizing Tests

You can modify the test parameters by editing `test_parallel_performance.R`:

```r
# Configuration (lines 14-15)
SEED <- 42                # Change random seed
TOLERANCE <- 1e-10        # Adjust numerical tolerance

# Test 1: Reproducibility (line 65)
test_reproducibility("MFA", mfa_em_fit, N = 50, K = 2, r = 2)
# Change N, K, r to test different configurations

# Test 2: Performance scaling (line 106)
test_performance_scaling("MFA", mfa_em_fit, N = 100, K = 2, r = 2)
# Change N to test with larger/smaller datasets

# Test 3: Data size scaling (line 148)
data_sizes <- c(20, 50, 100)
# Add more data sizes or change existing ones
```

## Integration with CI/CD

To integrate these tests into a CI/CD pipeline:

```bash
# Run performance tests
Rscript tests/performance/test_parallel_performance.R

# Check exit code
if [ $? -eq 0 ]; then
  echo "Performance tests passed"
else
  echo "Performance tests failed"
  exit 1
fi
```

## Troubleshooting

### Tests are slow

Performance tests involve actual timing measurements and may take 2-5 minutes to complete. This is expected.

### Inconsistent speedup results

System load, CPU frequency scaling, and other background processes can affect timing measurements. Run tests on a quiet system for most consistent results.

### OpenMP not available

If OpenMP is not available, the tests will still pass but speedup will be minimal (close to 1.0×). Check:

```r
# Check if OpenMP is enabled
Sys.getenv("OMP_NUM_THREADS")

# Set thread count explicitly
Sys.setenv(OMP_NUM_THREADS = "4")
```

### Numerical precision issues

If reproducibility tests fail with very small differences (< 1e-8), this may be due to floating-point arithmetic variations. Consider adjusting `TOLERANCE` if needed.

## Related Documentation

- **Benchmark Suite**: `benchmarks/README.md` - Comprehensive parallel efficiency benchmark
- **Main README**: `README.md` - Package overview and installation
- **OpenMP Implementation**: `src/mfaEM.cpp`, `src/mpcaEM.cpp` - C++ parallelization code

## Contact

For questions or issues:
- Open an issue on the GitHub repository
- Reference PR #31 (C++コードの最適化と並列処理改善)
