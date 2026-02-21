// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

/************************************************************
 * dmvnormPCA_rowwiseCpp
 ************************************************************/
// [[Rcpp::export]]
arma::vec dmvnormPCA_rowwiseCpp(const arma::mat &Xmat,
                                const arma::vec &mu,
                                const arma::mat &W,
                                double sigma2)
{
  int T = Xmat.n_rows;
  int M = Xmat.n_cols;

  // Sigma = W W^T + sigma2 I
  mat Sigma = W * W.t();
  Sigma.diag() += sigma2;
  // 数値安定
  Sigma.diag() += 1e-8;

  mat L = chol(Sigma, "lower");

  double logdetSigma = 0.0;
  for(int m=0; m<M; m++){
    logdetSigma += 2.0 * std::log(L(m,m));
  }
  double Mlog2pi = M * std::log(2.0 * M_PI);

  // invert L
  mat Linv = inv(trimatl(L));

  arma::vec out(T, fill::zeros);
  for(int i=0; i<T; i++){
    vec diff = Xmat.row(i).t() - mu;
    vec z = Linv * diff;
    double quadform = dot(z,z);
    double val = -0.5 * (Mlog2pi + logdetSigma + quadform);
    out[i] = val;
  }
  return out;
}

/************************************************************
 * pcaEMupdateCpp
 ************************************************************/
// [[Rcpp::export]]
Rcpp::List pcaEMupdateCpp(const arma::mat &Xmat, int r, int nIter)
{
  int n = Xmat.n_rows;
  int M = Xmat.n_cols;

  rowvec mu_r = mean(Xmat, 0);
  mat Xc = Xmat.each_row() - mu_r;

  // init W by SVD
  mat U, V;
  vec S;
  svd(U, S, V, Xc);
  mat Vr = V.cols(0, r-1);
  vec sr = S.subvec(0, r-1);
  mat W = Vr * diagmat(sr / std::sqrt((double)n));

  // init sigma^2
  double sigma2 = 1.0;
  if(r < M){
    double extraVar = 0.0;
    for(int j=r; j<M; j++){
      extraVar += std::pow(S[j], 2.0) / double(n);
    }
    extraVar /= double(M - r);
    sigma2 = extraVar;
  }

  // EM
  for(int it=0; it<nIter; it++){
    mat Sigma = W * W.t();
    Sigma.diag() += sigma2;
    Sigma.diag() += 1e-8;

    mat L = chol(Sigma, "lower");
    mat invL = inv(trimatl(L));
    mat invSigma = invL.t() * invL;  // (M x M)

    mat WtSi = W.t() * invSigma;     // (r x M)
    mat Mmat = eye(r, r) + WtSi * W; // (r x r)
    mat M_chol = chol(Mmat, "lower");
    mat M_inv = inv(trimatu(M_chol));
    M_inv = M_inv * M_inv.t();

    mat EZ = M_inv * WtSi * Xc.t();  // (r x n)
    mat EZZt = (EZ * EZ.t()) / double(n);
    EZZt += M_inv;

    mat XcEZt = (Xc.t() * EZ.t()) / double(n);
    mat EZZt_inv = inv(EZZt + 1e-12 * eye(r,r));
    mat W_new = XcEZt * EZZt_inv;

    // sigma2 update
    mat Sxx = (Xc.t() * Xc) / double(n);
    mat tmp = W_new * XcEZt.t();
    mat residual = Sxx - tmp;
    double traceVal = accu(residual.diag());
    double new_sigma2 = traceVal / double(M);
    if(new_sigma2 < 1e-12) new_sigma2 = 1e-12;

    W = W_new;
    sigma2 = new_sigma2;
  }

  Rcpp::List out;
  out["mu"]     = mu_r.t();
  out["W"]      = W;
  out["sigma2"] = sigma2;
  return out;
}

/************************************************************
 * pcaClosedFormCpp
 ************************************************************/
// [[Rcpp::export]]
Rcpp::List pcaClosedFormCpp(const arma::mat &Xmat, int r)
{
  int n = Xmat.n_rows;
  int M = Xmat.n_cols;

  rowvec mu_r = mean(Xmat, 0);
  mat Xc = Xmat.each_row() - mu_r;

  // SVD
  mat U, V;
  vec S;
  svd(U, S, V, Xc);

  mat Vr = V.cols(0, r-1);
  vec Sr = S.subvec(0, r-1);

  // leftover var => sigma2
  double sigma2 = 0.0;
  if(r < M){
    double sumLeft = 0.0;
    for(int j=r; j<M; j++){
      double lam_j = (S[j]*S[j]) / double(n);
      sumLeft += lam_j;
    }
    sigma2 = sumLeft / double(M - r);
  } else {
    sigma2 = 1e-12;
  }

  // W
  mat W = arma::zeros(M, r);
  for(int j=0; j<r; j++){
    double lam_j = (Sr[j]*Sr[j]) / double(n);
    double diff_j = lam_j - sigma2;
    if(diff_j < 1e-12) diff_j = 1e-12;
    double scale_j = std::sqrt(diff_j);
    W.col(j) = Vr.col(j) * scale_j;
  }

  Rcpp::List out;
  out["mu"]     = mu_r.t();
  out["W"]      = W;
  out["sigma2"] = sigma2;
  return out;
}

/************************************************************
 * pcaEMupdateWeightedCpp: Weighted PCA EM update
 * Used by soft EM for responsibility-weighted parameter estimation
 ************************************************************/
// [[Rcpp::export]]
Rcpp::List pcaEMupdateWeightedCpp(const arma::mat &Xmat, int r, int nIter,
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

  // Initial W via eigen decomposition of S
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
  arma::mat W = Vr * arma::diagmat(arma::sqrt(sr));

  // Initial sigma^2
  double sigma2 = 0.0;
  if(r < M){
    for(int j = r; j < M; j++){
      int jj = idx[j];
      sigma2 += std::max(eigval[jj], 0.0);
    }
    sigma2 /= double(M - r);
  }
  if(sigma2 < 1e-12) sigma2 = 1e-12;

  // EM iterations for weighted PCA
  for(int it = 0; it < nIter; it++){
    arma::mat Sigma = W * W.t();
    Sigma.diag() += sigma2;
    Sigma.diag() += 1e-8;

    // Compute Sigma inverse via Cholesky
    arma::mat L_chol = arma::chol(Sigma, "lower");
    arma::mat Linv = arma::inv(arma::trimatl(L_chol));
    arma::mat Siginv = Linv.t() * Linv;

    arma::mat WtSi = W.t() * Siginv;     // (r x M)
    arma::mat Mmat = arma::eye(r, r) + WtSi * W; // (r x r)
    arma::mat M_chol = arma::chol(Mmat, "lower");
    arma::mat M_inv = arma::inv(arma::trimatu(M_chol));
    M_inv = M_inv * M_inv.t();

    // E-step: compute expected sufficient statistics
    arma::mat sumWEzz = arma::zeros(r, r);   // sum of w_i * E[zz'|x_i]
    arma::mat sumWXEz = arma::zeros(M, r);   // sum of w_i * x_i * E[z|x_i]'

    for(int i = 0; i < n; i++){
      arma::vec xi = Xc.row(i).t();
      arma::vec Ez = M_inv * WtSi * xi;  // E[z|x]
      arma::mat Ezz = M_inv + Ez * Ez.t();  // E[zz'|x]

      sumWEzz += weights[i] * Ezz;
      sumWXEz += weights[i] * (xi * Ez.t());
    }

    // M-step: W = (sum w_i x_i E[z|x_i]') * (sum w_i E[zz'|x_i])^{-1}
    arma::mat W_new = sumWXEz * arma::inv(sumWEzz + 1e-12 * arma::eye(r, r));

    // Update sigma2: sigma2 = (1/M) * trace(S - W_new * sumWXEz' / wsum)
    arma::mat tmp = W_new * sumWXEz.t() / wsum;
    arma::mat residual = S - tmp;
    double traceVal = arma::accu(residual.diag());
    double new_sigma2 = traceVal / double(M);
    if(new_sigma2 < 1e-12) new_sigma2 = 1e-12;

    W = W_new;
    sigma2 = new_sigma2;
  }

  Rcpp::List out;
  out["mu"]     = mu_r.t();
  out["W"]      = W;
  out["sigma2"] = sigma2;
  return out;
}

/************************************************************
 * mpcaTimeseriesSoftCpp: Soft EM for Mixture PCA
 *   Uses responsibility-weighted parameter updates
 ************************************************************/
// [[Rcpp::export]]
Rcpp::List mpcaTimeseriesSoftCpp(Rcpp::List list_of_data,
                                  int K, int r,
                                  int max_iter,
                                  int nIterPCA,
                                  double tol,
                                  std::string method,
                                  IntegerVector z_init,
                                  bool verbose,
                                  int n_threads,
                                  int seed)
{
  int N = list_of_data.size();
  if(N < 1) stop("No data in list_of_data");

  NumericMatrix firstMat = list_of_data[0];
  int M = firstMat.ncol();

#ifdef _OPENMP
  if(n_threads > 0) {
    omp_set_num_threads(n_threads);
  }
  if(verbose) {
    int actual_threads = omp_get_max_threads();
    Rcpp::Rcout << "Soft MPCA EM using OpenMP with " << actual_threads << " threads\n";
  }
#endif

  // Initialize parameters from z_init
  std::vector<arma::vec> mu(K);
  std::vector<arma::mat> W(K);
  std::vector<double> sigma2(K);
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
      W[k] = 0.01 * arma::randn<arma::mat>(M, r);
      sigma2[k] = 0.1;
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

    Rcpp::List pcaRes;
    if(method == "EM"){
      pcaRes = pcaEMupdateCpp(X_k, r, nIterPCA);
    } else {
      pcaRes = pcaClosedFormCpp(X_k, r);
    }
    mu[k] = as<arma::vec>(pcaRes["mu"]);
    W[k] = as<arma::mat>(pcaRes["W"]);
    sigma2[k] = as<double>(pcaRes["sigma2"]);
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

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for(int i = 0; i < N; i++){
      for(int k = 0; k < K; k++){
        arma::vec lvec = dmvnormPCA_rowwiseCpp(X_list[i], mu[k], W[k], sigma2[k]);
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
          Rcpp::Rcout << "Soft MPCA EM converged at iter=" << iterEM
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

    // Update mu, W, sigma2 for each cluster
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

      // Use weighted PCA update
      Rcpp::List pcaRes = pcaEMupdateWeightedCpp(X_all, r, nIterPCA, w_all);
      mu[k] = as<arma::vec>(pcaRes["mu"]);
      W[k] = as<arma::mat>(pcaRes["W"]);
      sigma2[k] = as<double>(pcaRes["sigma2"]);
    }

    if(verbose && iterEM % 10 == 0){
      Rcpp::Rcout << "Soft MPCA EM iter " << iterEM << " loglik=" << new_loglik << "\n";
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
  Rcpp::List mu_list(K), W_list(K);
  NumericVector sigma2_vec(K), pi_out(K);
  for(int k = 0; k < K; k++){
    mu_list[k] = mu[k];
    W_list[k] = W[k];
    sigma2_vec[k] = sigma2[k];
    pi_out[k] = pi_k[k];
  }

  return Rcpp::List::create(
    Named("z") = z_final,
    Named("pi") = pi_out,
    Named("mu") = mu_list,
    Named("W") = W_list,
    Named("sigma2") = sigma2_vec,
    Named("gamma") = gamma
  );
}

/************************************************************
 * mpcaTimeseriesCpp:
 *   ハード割当Mixture PCA
 *   初期クラスタ z_init があれば使用、なければ (i%K)+1
 ************************************************************/
// [[Rcpp::export]]
Rcpp::List mpcaTimeseriesCpp(Rcpp::List list_of_data,
                             int K,
                             int r,
                             int max_iter,
                             int nIterPCA,
                             double tol,
                             std::string method = "EM",
                             Rcpp::Nullable<Rcpp::IntegerVector> z_init = R_NilValue,
                             bool verbose = true,
                             int n_threads = 0,
                             int seed = 0,
                             Rcpp::Nullable<Rcpp::Function> progress_callback = R_NilValue)
{
  int N = list_of_data.size();
  if(N < 1) stop("No data in list_of_data");

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

  // get M from the first item
  NumericMatrix firstMat = list_of_data[0];
  int M = firstMat.ncol();

  // initial cluster assignment
  IntegerVector zR(N);
  if(z_init.isNotNull()){
    Rcpp::IntegerVector zz(z_init.get());
    if(int(zz.size()) != N){
      stop("z_init size != N");
    }
    for(int i=0; i<N; i++){
      if(zz[i] < 1 || zz[i] > K){
        stop("z_init has out-of-range cluster label");
      }
      zR[i] = zz[i];
    }
  } else {
    for(int i=0; i<N; i++){
      zR[i] = (i % K) + 1;
    }
  }

  // parameters
  std::vector<arma::mat> W(K);
  std::vector<arma::vec> mu(K);
  std::vector<double> sigma2(K);

  // initialize using Armadillo RNG (seeded if seed != 0)
  if(seed != 0) {
    arma::arma_rng::set_seed(static_cast<unsigned int>(seed));
  }
  GetRNGstate();
  for(int k=0; k<K; k++){
    W[k] = 0.01 * arma::randn<arma::mat>(M, r);
    mu[k]= arma::vec(M, fill::zeros);
    sigma2[k] = 0.1;
  }
  PutRNGstate();

  arma::vec pi_k(K, fill::value(1.0 / double(K)));
  double old_loglik = -1e15;

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

    // ============ M-step ============
    for(int k2=0; k2<K; k2++){
      std::vector<int> idx;
      idx.reserve(N);
      for(int i=0; i<N; i++){
        if(zR[i] == (k2+1)) {
          idx.push_back(i);
        }
      }
      if(idx.empty()) {
        // skip
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

      // single-cluster PPCA
      Rcpp::List pcaRes;
      if(method == "EM"){
        pcaRes = pcaEMupdateCpp(X_k, r, nIterPCA);
      } else if(method=="closed_form"){
        pcaRes = pcaClosedFormCpp(X_k, r);
      } else {
        Rcpp::stop("method must be 'EM' or 'closed_form'.");
      }

      mu[k2]     = as<arma::vec>(pcaRes["mu"]);
      W[k2]      = as<arma::mat>(pcaRes["W"]);
      sigma2[k2] = as<double>(pcaRes["sigma2"]);
    }

    // ============ E-step ============
    arma::mat logLik_i(N, K, fill::zeros);
    arma::vec logPi(K, fill::zeros);
    for(int k2=0; k2<K; k2++){
      logPi[k2] = std::log(pi_k[k2] + 1e-16);
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
          arma::vec lvec = dmvnormPCA_rowwiseCpp(X_list[i], mu[k2], W[k2], sigma2[k2]);
          double sumLog = arma::accu(lvec);
          logLik_i(i, k2) = logPi[k2] + sumLog;
        }
      }
    }
#else
    for(int i=0; i<N; i++){
      for(int k2=0; k2<K; k2++){
        arma::vec lvec = dmvnormPCA_rowwiseCpp(X_list[i], mu[k2], W[k2], sigma2[k2]);
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
        double val = logLik_i(i, k2);
        if(val > bestVal){
          bestVal = val;
          bestk   = k2+1;
        }
      }
      zR_new[i] = bestk;
    }

    // pi_k update
    for(int k2=0; k2<K; k2++){
      int countk=0;
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

    // check for assignment convergence
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

  // final output
  Rcpp::IntegerVector z_out(N);
  for(int i=0; i<N; i++){
    z_out[i] = zR[i];
  }

  Rcpp::List W_out(K), Mu_out(K);
  NumericVector sigma2_out(K), pi_out(K);
  for(int k=0; k<K; k++){
    W_out[k]       = wrap(W[k]);
    Mu_out[k]      = wrap(mu[k]);
    sigma2_out[k]  = sigma2[k];
    pi_out[k]      = pi_k[k];
  }

  Rcpp::List ret;
  ret["z"]      = z_out;
  ret["W"]      = W_out;
  ret["mu"]     = Mu_out;
  ret["sigma2"] = sigma2_out;
  ret["pi"]     = pi_out;

  return ret;
}
