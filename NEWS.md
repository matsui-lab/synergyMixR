# synergyMixR 1.0.0

## Initial Public Release

* Core mixture model implementations: MPCA (Mixture of Probabilistic Component Analyzers) and MFA (Mixture of Factor Analyzers)
* EM algorithm implementations in C++ via Rcpp/RcppArmadillo for performance
* Soft EM with responsibility-weighted updates
* K-means based initialization for cluster assignments
* Comprehensive model selection tools:
  - BIC-based grid search over (K, r)
  - ICL (Integrated Complete-data Likelihood) with temperature parameter
  - Cross-Validation with 1-SE rule
  - Bootstrap stability selection for K
  - Two-stage selection for stability
* r-selection framework: global, clusterwise, and hybrid modes
* Rotation and flip utilities for synergy alignment
* Complete EMG preprocessing pipeline (DC removal, bandpass, notch, rectification, envelope, normalization)
* Visualization functions for muscle synergy analysis with unified ggplot2 theme
* Simulation utilities for generating mixture data with optional AR(1) temporal noise
* HPC parallel computation support via future/PBS
* Comprehensive testing infrastructure with testthat
* Quickstart vignette and package documentation
