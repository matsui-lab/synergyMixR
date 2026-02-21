# synergyMixR

[![R-universe](https://matsui-lab.r-universe.dev/badges/synergyMixR)](https://matsui-lab.r-universe.dev/synergyMixR)

An R package for mixture-based muscle synergy analysis (Mixture Factor Analysis and Mixture Principal Component Analysis).

## Installation

```r
# Install from R-universe (recommended, includes pre-built binaries)
install.packages("synergyMixR",
  repos = c("https://matsui-lab.r-universe.dev", "https://cloud.r-project.org"))

# Install from GitHub
devtools::install_github("matsui-lab/synergyMixR")

# Or install from CRAN (when available)
install.packages("synergyMixR")
```

## System Requirements

This package requires C++ compilation. Ensure you have:
- **Windows**: [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
- **macOS**: Xcode Command Line Tools (`xcode-select --install`)
- **Linux**: `build-essential` and `r-base-dev` packages

## Quick Start

```r
library(synergyMixR)

# Generate synthetic data
sim_data <- generate_synergy_data(
  N = 30, K = 2, r = 3, M = 8,
  T_each = 100, seed = 123
)

# Fit MFA model
fit_mfa <- mfa_em_fit(
  list_of_data = sim_data$list_of_data,
  K = 2, r = 3,
  n_init = 3,
  use_kmeans_init = TRUE
)

# Visualize synergy loadings
plot_cluster_synergy_loadings_mfa(fit_mfa, plot_type = "heatmap")
```

For a complete tutorial, see:

```r
vignette("quickstart", package = "synergyMixR")
```

## Key Features

- **Mixture Models**: MFA and MPCA with C++ backend (Rcpp/RcppArmadillo)
- **Model Selection**: BIC, ICL, cross-validation, bootstrap stability
- **EMG Preprocessing**: Complete pipeline from raw EMG to analysis-ready data
- **Visualization**: Unified ggplot2-based plotting functions
- **Simulation**: Synthetic data generation with configurable noise models
- **HPC Support**: Parallel computation via future/PBS backends

## Vignettes

| Vignette | Description |
|----------|-------------|
| `quickstart` | Quick start guide |
| `tutorial` | Comprehensive tutorial |
| `model_selection_guide` | Complete K,r selection workflow |
| `r_selection_framework` | r-selection methods |
| `metrics_reference` | Understanding metrics |
| `simulation_cluster_separation` | Cluster separation analysis |
| `simulation_single_vs_mixture` | Single vs mixture model comparison |

## Citation

If you use this package in your research, please cite:

```
Matsui, Y. (2026). synergyMixR: Mixture-Based Muscle Synergy Analysis.
R package version 1.0.0. https://github.com/matsui-lab/synergyMixR
```

## License

GPL-3
