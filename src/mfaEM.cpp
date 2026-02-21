// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

//----------------------------------------------------
// dmvnormFA_rowwiseCpp: (シングルスレッド) 多次元正規対数密度
//----------------------------------------------------

// [[Rcpp::export]]
arma::vec dmvnormFA_rowwiseCpp(const arma::mat &Xmat,
                               const arma::vec &mu,
                               const arma::mat &Lambda,
                               const arma::vec &Psi_diag)
{
  int T = Xmat.n_rows;
  int M = Xmat.n_cols;

  // Sigma = Lambda * Lambda.t() + diag(Psi_diag)
  mat Sigma = Lambda * Lambda.t();
  for(int m=0; m<M; m++){
    Sigma(m,m) += Psi_diag[m];
  }
  // Numerical stabilization
  Sigma.diag() += 1e-8;

  // Chol(Sigma)
  mat L = chol(Sigma, "lower");
  double logdetSigma = 0.0;
  for(int m=0; m<M; m++){
    logdetSigma += 2.0 * std::log(L(m,m));
  }
  double log2pi = std::log(2.0 * M_PI);
  double Mlog2pi = M * log2pi;

  // L^-1
  mat Linv = inv(trimatl(L));

  // Output vector
  arma::vec out(T, fill::zeros);

  // Row-wise loop
  for(int i = 0; i < T; i++){
    vec diff = Xmat.row(i).t() - mu;
    vec z = Linv * diff;
    double quadform = dot(z,z);
    double val = -0.5 * (Mlog2pi + logdetSigma + quadform);
    out[i] = val;
  }

  return out;
}

//----------------------------------------------------
// faEMupdateCpp: (単一クラスタ) FactorAnalyzer EM
//----------------------------------------------------

// [[Rcpp::export]]
Rcpp::List faEMupdateCpp(const arma::mat &Xmat, int r, int nIterFA)
{
  int n = Xmat.n_rows;
  int M = Xmat.n_cols;

  rowvec mu_r = mean(Xmat, 0);
  mat Xc = Xmat.each_row() - mu_r;

  // Initial Lambda => via PCA (SVD)
  mat U, V;
  vec S;
  svd(U, S, V, Xc);
  mat Vr = V.cols(0, r-1);
  vec sr = S.subvec(0, r-1);
  mat Lambda = Vr * diagmat(sr / std::sqrt((double)n));

  // Initial psi
  vec psi(M, fill::ones);
  double extraVar = 0.0;
  if(r < M){
    for(int j=r; j<M; j++){
      extraVar += S[j];
    }
    extraVar /= double(M - r) * double(n);
  } else {
    extraVar = 1.0;
  }
  psi *= extraVar;

  // EM iterations
  for(int it=0; it<nIterFA; it++){
    mat Sig = Lambda * Lambda.t();
    for(int m=0; m<M; m++){
      Sig(m,m) += psi[m];
    }
    // Numerical stabilization
    Sig.diag() += 1e-8;

    mat L = chol(Sig, "lower");
    mat invL = inv(trimatl(L));
    mat invSigma = invL.t() * invL;

    mat W = Lambda.t() * invSigma;
    mat A = eye(r, r) + W * Lambda;
    mat A_chol = chol(A, "lower");
    mat A_inv  = inv(trimatl(A_chol));
    mat Vz     = A_inv.t() * A_inv;

    mat WX = W * Xc.t();
    mat EZ = Vz * WX;

    mat EZZt = (EZ * EZ.t()) / double(n);
    for(int j=0; j<r; j++){
      EZZt(j,j) += Vz(j,j);
    }

    mat XcEZt = (Xc.t() * EZ.t()) / double(n);
    mat Lambda_new = XcEZt * inv(EZZt + 1e-12*eye(r,r));

    mat Sxx = (Xc.t() * Xc) / double(n);
    mat tmp = Lambda_new * XcEZt.t();
    mat residual = Sxx - tmp;
    for(int m=0; m<M; m++){
      double val = residual(m,m);
      if(val < 1e-12) val = 1e-12;
      psi[m] = val;
    }

    Lambda = Lambda_new;
  }

  // Return
  List out;
  out["mu"]     = mu_r.t();
  out["Lambda"] = Lambda;
  out["Psi"]    = diagmat(psi);
  return out;
}

//----------------------------------------------------
// faEMupdateWeightedCpp: Weighted Factor Analysis EM update
// Used by soft EM for responsibility-weighted parameter estimation
//----------------------------------------------------

// [[Rcpp::export]]
Rcpp::List faEMupdateWeightedCpp(const arma::mat &Xmat, int r, int nIterFA,
                                  const arma::vec &weights)
{
  int n = Xmat.n_rows;
  int M = Xmat.n_cols;

  // Weighted mean
  double wsum = arma::accu(weights);
  arma::rowvec mu_r = arma::zeros<arma::rowvec>(M);
  for(int i = 0; i < n; i++){
    mu_r += weights[i] * Xmat.row(i);
  }
  mu_r /= wsum;

  // Weighted centered data
  arma::mat Xc = Xmat.each_row() - mu_r;

  // Weighted covariance
  arma::mat S = arma::zeros(M, M);
  for(int i = 0; i < n; i++){
    S += weights[i] * (Xc.row(i).t() * Xc.row(i));
  }
  S /= wsum;

  // Initial Lambda via eigen decomposition of S
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, S);

  // Sort eigenvalues descending
  arma::uvec idx = arma::sort_index(eigval, "descend");
  arma::mat Vr = eigvec.cols(idx.head(r));
  arma::vec sr = eigval(idx.head(r));
  for(int j = 0; j < r; j++){
    if(sr[j] < 1e-6) sr[j] = 1e-6;
  }
  arma::mat Lambda = Vr * arma::diagmat(arma::sqrt(sr));

  // Initial psi
  arma::vec psi(M, arma::fill::ones);
  double extraVar = 0.0;
  if(r < M){
    for(int j = r; j < M; j++){
      int jj = idx[j];
      extraVar += std::max(eigval[jj], 0.0);
    }
    extraVar /= double(M - r);
  } else {
    extraVar = 0.01;
  }
  psi *= extraVar;
  for(int m = 0; m < M; m++){
    if(psi[m] < 1e-6) psi[m] = 1e-6;
  }

  // EM iterations for weighted FA
  for(int it = 0; it < nIterFA; it++){
    arma::mat Sig = Lambda * Lambda.t();
    for(int m = 0; m < M; m++){
      Sig(m,m) += psi[m];
    }
    Sig.diag() += 1e-8;

    // Compute Sigma inverse via Cholesky
    arma::mat L_chol = arma::chol(Sig, "lower");
    arma::mat Linv = arma::inv(arma::trimatl(L_chol));
    arma::mat Siginv = Linv.t() * Linv;

    arma::mat Beta = Lambda.t() * Siginv;  // r x M

    // E-step: compute expected sufficient statistics
    arma::mat I_r = arma::eye(r, r);
    arma::mat Vz = I_r - Beta * Lambda;  // E[zz'|x] - E[z|x]E[z|x]' (constant part)

    arma::mat sumWEzz = arma::zeros(r, r);   // sum of w_i * E[zz'|x_i]
    arma::mat sumWEzX = arma::zeros(r, M);   // sum of w_i * E[z|x_i] * x_i'

    for(int i = 0; i < n; i++){
      arma::vec xi = Xc.row(i).t();
      arma::vec Ez = Beta * xi;  // E[z|x]
      arma::mat Ezz = Vz + Ez * Ez.t();  // E[zz'|x]

      sumWEzz += weights[i] * Ezz;
      sumWEzX += weights[i] * (Ez * xi.t());
    }

    // M-step: Lambda = (sum w_i x_i E[z|x_i]') * (sum w_i E[zz'|x_i])^{-1}
    arma::mat sumWXEz = sumWEzX.t();  // M x r
    arma::mat Lambda_new = sumWXEz * arma::inv(sumWEzz + 1e-12 * I_r);

    // Update psi: psi_m = (1/wsum) * sum w_i (x_im - mu_m)^2 - Lambda_new[m,:] * E[zx']_m / wsum
    // Simplified: use residual variance
    arma::mat Sxx_weighted = S;  // Already computed weighted covariance
    arma::mat tmp = Lambda_new * sumWEzX / wsum;
    for(int m = 0; m < M; m++){
      double val = Sxx_weighted(m,m) - tmp(m,m);
      if(val < 1e-6) val = 1e-6;
      psi[m] = val;
    }

    Lambda = Lambda_new;
  }

  arma::mat Psi = arma::diagmat(psi);
  arma::vec mu_out = mu_r.t();

  return Rcpp::List::create(
    Named("mu") = mu_out,
    Named("Lambda") = Lambda,
    Named("Psi") = Psi
  );
}

//----------------------------------------------------
// mfaTimeseriesSoftCpp: Soft EM for Mixture of Factor Analyzers
//    - Uses responsibility-weighted parameter updates
//    - Returns both hard assignments and soft responsibilities
//----------------------------------------------------

// [[Rcpp::export]]
Rcpp::List mfaTimeseriesSoftCpp(Rcpp::List list_of_data,
                                 int K, int r,
                                 int max_iter,
                                 int nIterFA,
                                 double tol,
                                 IntegerVector z_init,
                                 bool verbose,
                                 int n_threads,
                                 int seed)
{
  int N = list_of_data.size();
  if(N < 1) {
    stop("No data in list_of_data");
  }

  NumericMatrix firstMat = list_of_data[0];
  int M = firstMat.ncol();

#ifdef _OPENMP
  if(n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
  if(verbose) {
    int actual_threads = omp_get_max_threads();
    Rcpp::Rcout << "Soft EM using OpenMP with " << actual_threads << " threads\n";
  }
#endif

  // Initialize parameters from z_init
  std::vector<arma::vec> mu(K);
  std::vector<arma::mat> Lambda(K);
  std::vector<arma::mat> Psi(K);
  arma::vec pi_k(K, arma::fill::zeros);

  // Count initial assignments
  for(int i = 0; i < N; i++){
    pi_k[z_init[i] - 1] += 1.0;
  }
  pi_k /= double(N);

  // Initialize cluster parameters from initial assignment
  for(int k = 0; k < K; k++){
    std::vector<int> idx;
    for(int i = 0; i < N; i++){
      if(z_init[i] == k + 1) idx.push_back(i);
    }
    if(idx.empty()){
      // Empty cluster - initialize with random small values
      mu[k] = arma::zeros<arma::vec>(M);
      Lambda[k] = 0.01 * arma::randn<arma::mat>(M, r);
      Psi[k] = arma::eye<arma::mat>(M, M) * 0.1;
      continue;
    }

    long totalT = 0;
    for(size_t ii = 0; ii < idx.size(); ii++){
      NumericMatrix mat_i = list_of_data[idx[ii]];
      totalT += mat_i.nrow();
    }
    arma::mat X_k(totalT, M);
    long rowpos = 0;
    for(size_t ii = 0; ii < idx.size(); ii++){
      NumericMatrix mat_i = list_of_data[idx[ii]];
      int Ti = mat_i.nrow();
      arma::mat tmp(mat_i.begin(), Ti, M, false);
      X_k.rows(rowpos, rowpos + Ti - 1) = tmp;
      rowpos += Ti;
    }

    Rcpp::List fares = faEMupdateCpp(X_k, r, nIterFA);
    mu[k] = as<arma::vec>(fares["mu"]);
    Lambda[k] = as<arma::mat>(fares["Lambda"]);
    Psi[k] = as<arma::mat>(fares["Psi"]);
  }

  // Responsibility matrix
  arma::mat gamma(N, K, arma::fill::zeros);
  double old_loglik = -1e15;

  // Pre-convert data to arma::mat for efficiency
  std::vector<arma::mat> X_list(N);
  for(int i = 0; i < N; i++){
    NumericMatrix mat_i = list_of_data[i];
    X_list[i] = Rcpp::as<arma::mat>(mat_i);
  }

  // EM iterations
  for(int iterEM = 0; iterEM < max_iter; iterEM++){

    // E-step: compute responsibilities
    arma::mat logLik_i(N, K, arma::fill::zeros);
    arma::vec logPi(K);
    for(int k = 0; k < K; k++){
      logPi[k] = std::log(pi_k[k] + 1e-16);
    }

    std::vector<arma::vec> psi_diag_vec(K);
    for(int k = 0; k < K; k++){
      psi_diag_vec[k] = arma::vec(M);
      for(int m = 0; m < M; m++){
        psi_diag_vec[k][m] = Psi[k](m, m);
      }
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for(int i = 0; i < N; i++){
      for(int k = 0; k < K; k++){
        arma::vec lvec = dmvnormFA_rowwiseCpp(X_list[i], mu[k], Lambda[k], psi_diag_vec[k]);
        double sumLog = arma::accu(lvec);
        logLik_i(i, k) = logPi[k] + sumLog;
      }
    }

    // Compute responsibilities with log-sum-exp
    double new_loglik = 0.0;
    for(int i = 0; i < N; i++){
      double maxLog = logLik_i.row(i).max();
      double sumExp = 0.0;
      for(int k = 0; k < K; k++){
        sumExp += std::exp(logLik_i(i, k) - maxLog);
      }
      double logSumExp = maxLog + std::log(sumExp);
      new_loglik += logSumExp;

      for(int k = 0; k < K; k++){
        gamma(i, k) = std::exp(logLik_i(i, k) - logSumExp);
      }
    }

    // Check convergence
    if(iterEM >= 5){
      double diff = std::fabs(new_loglik - old_loglik);
      if(diff < tol){
        if(verbose){
          Rcpp::Rcout << "Soft EM converged at iter=" << iterEM
                      << " diff=" << diff << "\n";
        }
        break;
      }
    }
    old_loglik = new_loglik;

    // M-step: update parameters with responsibilities

    // Update pi
    for(int k = 0; k < K; k++){
      pi_k[k] = arma::accu(gamma.col(k)) / double(N);
    }

    // Update mu, Lambda, Psi for each cluster
    for(int k = 0; k < K; k++){
      arma::vec weights_k = gamma.col(k);
      double Nk = arma::accu(weights_k);

      if(Nk < 1e-10) continue;  // Skip empty cluster

      // Concatenate all time series data with observation-level weights
      long totalT = 0;
      for(int i = 0; i < N; i++){
        totalT += X_list[i].n_rows;
      }

      arma::mat X_all(totalT, M);
      arma::vec w_all(totalT);
      long pos = 0;
      for(int i = 0; i < N; i++){
        int Ti = X_list[i].n_rows;
        X_all.rows(pos, pos + Ti - 1) = X_list[i];
        for(int t = 0; t < Ti; t++){
          w_all[pos + t] = weights_k[i];
        }
        pos += Ti;
      }

      // Use weighted FA update
      Rcpp::List fares = faEMupdateWeightedCpp(X_all, r, nIterFA, w_all);
      mu[k] = as<arma::vec>(fares["mu"]);
      Lambda[k] = as<arma::mat>(fares["Lambda"]);
      Psi[k] = as<arma::mat>(fares["Psi"]);
    }

    if(verbose && iterEM % 10 == 0){
      Rcpp::Rcout << "Soft EM iter " << iterEM << " loglik=" << new_loglik << "\n";
    }

    // Check for user interrupt
    if(iterEM % 5 == 0){
      Rcpp::checkUserInterrupt();
    }
  }

  // Hard assignment from final responsibilities
  IntegerVector z_final(N);
  for(int i = 0; i < N; i++){
    int bestk = 0;
    double bestval = gamma(i, 0);
    for(int k = 1; k < K; k++){
      if(gamma(i, k) > bestval){
        bestval = gamma(i, k);
        bestk = k;
      }
    }
    z_final[i] = bestk + 1;
  }

  // Convert to R list format
  Rcpp::List mu_list(K), Lambda_list(K), Psi_list(K);
  for(int k = 0; k < K; k++){
    mu_list[k] = mu[k];
    Lambda_list[k] = Lambda[k];
    Psi_list[k] = Psi[k];
  }

  NumericVector pi_out(K);
  for(int k = 0; k < K; k++){
    pi_out[k] = pi_k[k];
  }

  return Rcpp::List::create(
    Named("z") = z_final,
    Named("pi") = pi_out,
    Named("mu") = mu_list,
    Named("Lambda") = Lambda_list,
    Named("Psi") = Psi_list,
    Named("gamma") = gamma
  );
}

//----------------------------------------------------
// mfaTimeseriesCpp: Mixture of Factor Analyzers (hard assignment)
//    - If z_init is provided, use it as the initial assignment
//    - Otherwise default to (i % K) + 1
//----------------------------------------------------

// [[Rcpp::export]]
Rcpp::List mfaTimeseriesCpp(Rcpp::List list_of_data,
                            int K,
                            int r,
                            int max_iter,
                            int nIterFA,
                            double tol,
                            Rcpp::Nullable<Rcpp::IntegerVector> z_init = R_NilValue,
                            bool verbose = true,
                            int n_threads = 0,
                            int seed = 0,
                            Rcpp::Nullable<Rcpp::Function> progress_callback = R_NilValue)
{
  int N = list_of_data.size();
  if(N < 1) {
    stop("No data in list_of_data");
  }

#ifdef _OPENMP
  if(n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
  int actual_threads = omp_get_max_threads();
  if(verbose) {
    if(seed != 0) {
      Rcpp::Rcout << "Using OpenMP with " << actual_threads << " threads (seed=" << seed << ")\n";
    } else {
      Rcpp::Rcout << "Using OpenMP with " << actual_threads << " threads\n";
    }
  }
#else
  if(verbose && n_threads > 0) {
    Rcpp::Rcout << "OpenMP not available, running single-threaded\n";
  }
#endif

  // Get M from the first subject
  NumericMatrix firstMat = list_of_data[0];
  int M = firstMat.ncol();

  // Initialize cluster labels zR
  IntegerVector zR(N);

  // If user provided an initial assignment
  if(z_init.isNotNull()){
    Rcpp::IntegerVector zz(z_init.get());
    if(zz.size() != N){
      stop("z_init size != N");
    }
    // Check 1..K
    for(int i=0; i<N; i++){
      if(zz[i] < 1 || zz[i] > K){
        stop("z_init has out-of-range cluster label");
      }
      zR[i] = zz[i];
    }
  } else {
    // Default init
    for(int i=0; i<N; i++){
      zR[i] = (i % K) + 1;
    }
  }

  // Parameters
  std::vector<arma::mat> Lambda(K);
  std::vector<arma::vec> mu(K);
  std::vector<arma::mat> Psi(K);

  // Random init using Armadillo RNG (seeded if seed != 0)
  if(seed != 0) {
    arma::arma_rng::set_seed(static_cast<unsigned int>(seed));
  }
  GetRNGstate();
  for(int k=0; k<K; k++){
    Lambda[k] = 0.01 * arma::randn<arma::mat>(M, r);
    mu[k]     = arma::vec(M, fill::zeros);
    Psi[k]    = arma::eye<arma::mat>(M,M)*0.1;
  }
  PutRNGstate();

  arma::vec pi_k(K, fill::value(1.0 / double(K)));
  double old_loglik = -1e15;

  // EM loop
  for(int iterEM = 1; iterEM <= max_iter; iterEM++){
    if(progress_callback.isNotNull()) {
#ifdef _OPENMP
#pragma omp master
#endif
      {
        Rcpp::Function callback(progress_callback);
        callback(iterEM, max_iter);
      }
    }
    
    if(verbose && iterEM % 10 == 0) {
      Rcpp::Rcout << "EM iteration " << iterEM << "/" << max_iter << "\n";
      Rcpp::checkUserInterrupt();
    }
    //==============================
    // M-step
    //==============================
    for(int k=0; k<K; k++){
      // gather indices for cluster k
      std::vector<int> idx;
      idx.reserve(N);
      for(int i=0; i<N; i++){
        if(zR[i] == (k+1)){
          idx.push_back(i);
        }
      }
      if(idx.empty()) {
        // skip empty cluster
        continue;
      }

      long totalT = 0;
      for(size_t ii=0; ii<idx.size(); ii++){
        NumericMatrix mat_i = list_of_data[idx[ii]];
        totalT += mat_i.nrow();
      }
      arma::mat X_k(totalT, M);
      long rowpos = 0;
      for(size_t ii=0; ii<idx.size(); ii++){
        NumericMatrix mat_i = list_of_data[idx[ii]];
        int Ti = mat_i.nrow();
        arma::mat tmp(mat_i.begin(), Ti, M, false);
        X_k.rows(rowpos, rowpos + Ti - 1) = tmp;
        rowpos += Ti;
      }

      // single-cluster Factor Analysis
      Rcpp::List fares = faEMupdateCpp(X_k, r, nIterFA);
      mu[k]     = as<arma::vec>(fares["mu"]);
      Lambda[k] = as<arma::mat>(fares["Lambda"]);
      Psi[k]    = as<arma::mat>(fares["Psi"]);
    }

    //==============================
    // E-step (hard assignment)
    //==============================
    arma::mat logLik_i(N, K, fill::zeros);
    arma::vec logPi(K);
    for(int k2=0; k2<K; k2++){
      logPi[k2] = std::log(pi_k[k2] + 1e-16);
    }

    std::vector<arma::vec> psi_diag_vec(K);
    for(int k2=0; k2<K; k2++){
      psi_diag_vec[k2] = arma::vec(M);
      for(int mm=0; mm<M; mm++){
        psi_diag_vec[k2][mm] = Psi[k2](mm,mm);
      }
    }

    std::vector<arma::mat> X_list(N);
    for(int i=0; i<N; i++){
      NumericMatrix mat_i = list_of_data[i];
      X_list[i] = Rcpp::as<arma::mat>(mat_i);
    }

#ifdef _OPENMP
#pragma omp parallel
    {
      if(seed != 0) {
        arma::arma_rng::set_seed(static_cast<unsigned int>(seed) + omp_get_thread_num());
      }
#pragma omp for schedule(dynamic)
      for(int i=0; i<N; i++){
        for(int k2=0; k2<K; k2++){
          arma::vec lvec = dmvnormFA_rowwiseCpp(X_list[i], mu[k2], Lambda[k2], psi_diag_vec[k2]);
          double sumLog = arma::accu(lvec);
          logLik_i(i, k2) = logPi[k2] + sumLog;
        }
      }
    }
#else
    for(int i=0; i<N; i++){
      for(int k2=0; k2<K; k2++){
        arma::vec lvec = dmvnormFA_rowwiseCpp(X_list[i], mu[k2], Lambda[k2], psi_diag_vec[k2]);
        double sumLog = arma::accu(lvec);
        logLik_i(i, k2) = logPi[k2] + sumLog;
      }
    }
#endif

    // Hard assignment
    arma::vec rowMax = arma::max(logLik_i, 1);
    IntegerVector zR_new(N);
    for(int i=0; i<N; i++){
      double bestVal = -1e15;
      int bestk = 1;
      for(int k2=0; k2<K; k2++){
        double val = logLik_i(i,k2);
        if(val > bestVal){
          bestVal = val;
          bestk   = k2 + 1;
        }
      }
      zR_new[i] = bestk;
    }

    // update pi_k
    for(int k2=0; k2<K; k2++){
      int countk = 0;
      for(int i=0; i<N; i++){
        if(zR_new[i] == (k2+1)){
          countk++;
        }
      }
      pi_k[k2] = double(countk) / double(N);
    }

    // log-likelihood
    double new_loglik = arma::accu(rowMax);
    double diff = std::fabs(new_loglik - old_loglik);

    // Check convergence
    // 1) If assignments haven't changed, stop
    bool sameAll = true;
    for(int i=0; i<N; i++){
      if(zR[i] != zR_new[i]){
        sameAll = false;
        break;
      }
    }

    if(iterEM >= 5){
      if(sameAll){
        Rcpp::Rcout << "Converged by assignment at iter=" << iterEM << "\n";
        zR = zR_new;
        break;
      }

      // 2) If log-likelihood diff < tol, stop
      if(diff < tol){
        Rcpp::Rcout << "Converged by loglik diff at iter=" << iterEM
                    << " diff=" << diff << "\n";
        zR = zR_new;
        break;
      }
    }
    zR = zR_new;
    old_loglik = new_loglik;
  }

  // Prepare output
  Rcpp::IntegerVector z_out(N);
  for(int i=0; i<N; i++){
    z_out[i] = zR[i];
  }

  Rcpp::List Lmbd(K), Mu(K), PsiList(K);
  for(int k=0; k<K; k++){
    Lmbd[k]    = wrap(Lambda[k]);
    Mu[k]      = wrap(mu[k]);
    PsiList[k] = wrap(Psi[k]);
  }

  NumericVector pi_out(K);
  for(int k2=0; k2<K; k2++){
    pi_out[k2] = pi_k[k2];
  }

  Rcpp::List ret;
  ret["z"]      = z_out;
  ret["Lambda"] = Lmbd;
  ret["mu"]     = Mu;
  ret["Psi"]    = PsiList;
  ret["pi"]     = pi_out;

  return ret;
}
