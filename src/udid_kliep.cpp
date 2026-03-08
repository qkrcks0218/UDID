// udid_kliep.cpp
//
// Fast replacement for compute_density_ratio_Kernel() in the UDiD estimator.
//
// Replaces the nested apply() R loops in densratio:::compute_kernel_Gaussian
// with a tight C++ double loop.
//
// OpenMP parallelism is used automatically when available (Linux/GCC).
// On macOS with Apple clang, it compiles and runs correctly on a single
// thread without any extra configuration. To enable OpenMP on macOS:
//   1. brew install libomp
//   2. Create ~/.R/Makevars containing:
//        CXXFLAGS += -Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include
//        LDFLAGS  += -L$(brew --prefix libomp)/lib -lomp
//
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cmath>

using namespace Rcpp;
using namespace arma;

//' Fast Gaussian kernel density ratio evaluation
//'
//' Computes  r(x_i) = sum_c w_c * exp(-||x_i - center_c||^2 / (2*sigma^2))
//' for every row x_i of \code{x}.
//'
//' Drop-in replacement for:
//'   as.vector(densratio:::compute_kernel_Gaussian(x, centers, sigma) %*% kernel_weights)
//'
//' @param x              (N_eval x p)    evaluation points
//' @param centers        (N_centers x p) kernel centres from KLIEP fit
//' @param sigma          scalar          bandwidth
//' @param kernel_weights (N_centers)     weights from KLIEP fit
//' @param n_threads      int             OpenMP threads (0 = all; ignored without OpenMP)
//' @return numeric vector of length N_eval
// [[Rcpp::export]]
arma::vec gaussian_kernel_ratio(
    const arma::mat& x,
    const arma::mat& centers,
    double           sigma,
    const arma::vec& kernel_weights,
    int              n_threads = 0)
{
  const int    n_eval    = (int)x.n_rows;
  const int    n_centers = (int)centers.n_rows;
  const int    n_cols    = (int)x.n_cols;
  const double inv_2s2   = 1.0 / (2.0 * sigma * sigma);

  arma::vec result(n_eval, arma::fill::zeros);

#ifdef _OPENMP
  if (n_threads > 0) omp_set_num_threads(n_threads);
  #pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < n_eval; i++) {
    double r = 0.0;
    for (int c = 0; c < n_centers; c++) {
      double dist2 = 0.0;
      for (int p = 0; p < n_cols; p++) {
        double d = x(i, p) - centers(c, p);
        dist2 += d * d;
      }
      r += kernel_weights[c] * std::exp(-dist2 * inv_2s2);
    }
    result[i] = r;
  }

  return result;
}


//' Fast Gaussian kernel matrix with per-dimension bandwidths
//'
//' Computes  K_{i,c} = exp( - sum_j (x_{i,j} - center_{c,j})^2 / (2 * sigma_j^2) )
//' for every row i and every centre c.
//'
//' @param x        (N x p) evaluation points
//' @param centers  (M x p) kernel centres
//' @param sigma_vec (p)    per-dimension bandwidths
//' @param n_threads int    OpenMP threads (0 = all; ignored without OpenMP)
//' @return (N x M) kernel matrix
// [[Rcpp::export]]
arma::mat gaussian_kernel_matrix_bw(
    const arma::mat& x,
    const arma::mat& centers,
    const arma::vec& sigma_vec,
    int              n_threads = 0)
{
  const int n_eval    = (int)x.n_rows;
  const int n_centers = (int)centers.n_rows;
  const int n_cols    = (int)x.n_cols;

  // precompute 1/(2*sigma_j^2) for each dimension
  arma::vec inv_2s2(n_cols);
  for (int p = 0; p < n_cols; p++) {
    inv_2s2[p] = 1.0 / (2.0 * sigma_vec[p] * sigma_vec[p]);
  }

  arma::mat result(n_eval, n_centers);

#ifdef _OPENMP
  if (n_threads > 0) omp_set_num_threads(n_threads);
  #pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < n_eval; i++) {
    for (int c = 0; c < n_centers; c++) {
      double dist2 = 0.0;
      for (int p = 0; p < n_cols; p++) {
        double d = x(i, p) - centers(c, p);
        dist2 += d * d * inv_2s2[p];
      }
      result(i, c) = std::exp(-dist2);
    }
  }

  return result;
}


//' Fast Gaussian kernel density ratio evaluation with per-dimension bandwidths
//'
//' Computes  r(x_i) = sum_c w_c * exp( - sum_j (x_{i,j} - center_{c,j})^2 / (2 sigma_j^2) )
//'
//' @param x              (N_eval x p)    evaluation points
//' @param centers        (N_centers x p) kernel centres
//' @param sigma_vec      (p)             per-dimension bandwidths
//' @param kernel_weights (N_centers)     weights from KLIEP fit
//' @param n_threads      int             OpenMP threads (0 = all; ignored without OpenMP)
//' @return numeric vector of length N_eval
// [[Rcpp::export]]
arma::vec gaussian_kernel_ratio_bw(
    const arma::mat& x,
    const arma::mat& centers,
    const arma::vec& sigma_vec,
    const arma::vec& kernel_weights,
    int              n_threads = 0)
{
  const int n_eval    = (int)x.n_rows;
  const int n_centers = (int)centers.n_rows;
  const int n_cols    = (int)x.n_cols;

  arma::vec inv_2s2(n_cols);
  for (int p = 0; p < n_cols; p++) {
    inv_2s2[p] = 1.0 / (2.0 * sigma_vec[p] * sigma_vec[p]);
  }

  arma::vec result(n_eval, arma::fill::zeros);

#ifdef _OPENMP
  if (n_threads > 0) omp_set_num_threads(n_threads);
  #pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < n_eval; i++) {
    double r = 0.0;
    for (int c = 0; c < n_centers; c++) {
      double dist2 = 0.0;
      for (int p = 0; p < n_cols; p++) {
        double d = x(i, p) - centers(c, p);
        dist2 += d * d * inv_2s2[p];
      }
      r += kernel_weights[c] * std::exp(-dist2);
    }
    result[i] = r;
  }

  return result;
}


// ================================================================
//  Full KLIEP implementation in C++ for speed
// ================================================================

// Compute the kernel matrix (N x M) with per-dimension bandwidths
static arma::mat compute_phi(const arma::mat& x,
                              const arma::mat& centers,
                              const arma::vec& inv_2s2)
{
  const int N = (int)x.n_rows;
  const int M = (int)centers.n_rows;
  const int p = (int)x.n_cols;
  arma::mat phi(N, M);

  for (int i = 0; i < N; i++) {
    for (int c = 0; c < M; c++) {
      double d2 = 0.0;
      for (int j = 0; j < p; j++) {
        double d = x(i, j) - centers(c, j);
        d2 += d * d * inv_2s2[j];
      }
      phi(i, c) = std::exp(-d2);
    }
  }
  return phi;
}

// KLIEP score: mean(log(phi * alpha)), returns -1e30 on failure
static double kliep_score(const arma::mat& phi, const arma::vec& alpha)
{
  arma::vec vals = phi * alpha;
  double s = 0.0;
  int n = (int)vals.n_elem;
  for (int i = 0; i < n; i++) {
    if (vals[i] <= 0.0) return -1e30;
    s += std::log(vals[i]);
  }
  return s / n;
}

// Project alpha: alpha += (1 - b'alpha)*c, clip negatives, normalise by b'alpha
static void project_alpha(arma::vec& alpha, const arma::vec& b, const arma::vec& c_vec)
{
  double ba = arma::dot(b, alpha);
  alpha += (1.0 - ba) * c_vec;
  // clip negatives
  for (arma::uword i = 0; i < alpha.n_elem; i++) {
    if (alpha[i] < 0.0) alpha[i] = 0.0;
  }
  ba = arma::dot(b, alpha);
  if (ba > 0.0) alpha /= ba;
}

// KLIEP optimise alpha (matches densratio:::KLIEP_optimize_alpha)
static arma::vec kliep_optimize(const arma::mat& phi_x1,
                                 const arma::vec& mean_phi_x2)
{
  const int M = (int)phi_x1.n_cols;
  arma::vec b = mean_phi_x2;
  double bb = arma::dot(b, b);
  if (bb < 1e-30) return arma::vec(M, arma::fill::zeros);
  arma::vec c_vec = b / bb;

  arma::vec alpha(M, arma::fill::ones);
  project_alpha(alpha, b, c_vec);
  double score = kliep_score(phi_x1, alpha);

  double epsilons[] = {1000, 100, 10, 1, 0.1, 0.01, 0.001};
  const int max_iter = 100;

  for (int e = 0; e < 7; e++) {
    double eps = epsilons[e];
    for (int it = 0; it < max_iter; it++) {
      // gradient step:  alpha_new = alpha + eps * A^T * (1 / (A * alpha))
      arma::vec Aa = phi_x1 * alpha;
      arma::vec inv_Aa(Aa.n_elem);
      for (arma::uword i = 0; i < Aa.n_elem; i++) {
        inv_Aa[i] = (Aa[i] > 1e-30) ? 1.0 / Aa[i] : 1e30;
      }
      arma::vec alpha_new = alpha + eps * (phi_x1.t() * inv_Aa);
      project_alpha(alpha_new, b, c_vec);
      double score_new = kliep_score(phi_x1, alpha_new);
      if (score_new <= score) break;
      alpha = alpha_new;
      score = score_new;
    }
  }
  return alpha;
}

// CV score for a given sigma_vec: build kernel matrices, run K-fold KLIEP
static double kliep_cv_score(const arma::mat& x1,
                              const arma::mat& x2,
                              const arma::mat& centers,
                              const arma::vec& inv_2s2,
                              const arma::ivec& cv_split,
                              int fold)
{
  arma::mat phi_x1 = compute_phi(x1, centers, inv_2s2);
  arma::mat phi_x2 = compute_phi(x2, centers, inv_2s2);

  // check for degenerate kernels
  if (!phi_x1.is_finite() || !phi_x2.is_finite()) return -1e30;

  arma::vec mean_phi_x2 = arma::mean(phi_x2, 0).t();  // column means
  if (arma::accu(arma::abs(mean_phi_x2)) < 1e-30) return -1e30;

  double total = 0.0;
  for (int k = 1; k <= fold; k++) {
    // extract train/test indices
    arma::uvec train_idx = arma::find(cv_split != k);
    arma::uvec test_idx  = arma::find(cv_split == k);
    if (train_idx.n_elem == 0 || test_idx.n_elem == 0) return -1e30;

    arma::vec alpha = kliep_optimize(phi_x1.rows(train_idx), mean_phi_x2);
    double s = kliep_score(phi_x1.rows(test_idx), alpha);
    if (s <= -1e29) return -1e30;
    total += s;
  }
  return total / fold;
}

// Single-coordinate search: mirrors .search_one() in R
static double search_one_coord(double sigma_init,
                                const arma::mat& x1,
                                const arma::mat& x2,
                                const arma::mat& centers,
                                arma::vec& sigma_vec,  // modified in-place
                                const arma::uvec& coord_idx,
                                const arma::ivec& cv_split,
                                int fold)
{
  auto set_and_score = [&](double s) -> double {
    arma::vec inv_2s2(sigma_vec.n_elem);
    for (arma::uword j = 0; j < sigma_vec.n_elem; j++) {
      double sv = sigma_vec[j];
      inv_2s2[j] = 1.0 / (2.0 * sv * sv);
    }
    // temporarily set the coordinates
    for (arma::uword ci = 0; ci < coord_idx.n_elem; ci++) {
      double tmp = 1.0 / (2.0 * s * s);
      inv_2s2[coord_idx[ci]] = tmp;
    }
    return kliep_cv_score(x1, x2, centers, inv_2s2, cv_split, fold);
  };

  double sigma = sigma_init;
  double score = set_and_score(sigma);
  if (!std::isfinite(score)) score = -1e30;

  double digit_powers[] = {10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001};
  for (int dp = 0; dp < 7; dp++) {
    double step = digit_powers[dp];
    // try increasing
    for (int i = 0; i < 9; i++) {
      double sigma_new = sigma + step;
      if (sigma_new <= 0) break;
      double score_new = set_and_score(sigma_new);
      if (!std::isfinite(score_new) || score_new <= score) break;
      score = score_new;
      sigma = sigma_new;
    }
    // try decreasing
    for (int i = 0; i < 9; i++) {
      double sigma_new = sigma - step;
      if (sigma_new <= 1e-6) break;
      double score_new = set_and_score(sigma_new);
      if (!std::isfinite(score_new) || score_new <= score) break;
      score = score_new;
      sigma = sigma_new;
    }
  }

  // update sigma_vec with the found value
  for (arma::uword ci = 0; ci < coord_idx.n_elem; ci++) {
    sigma_vec[coord_idx[ci]] = sigma;
  }
  return sigma;
}


//' Full KLIEP with per-dimension bandwidth selection (all in C++)
//'
//' Runs coordinate-descent sigma search + KLIEP weight optimisation.
//' Drop-in replacement for KLIEP_bw() in R but much faster.
//'
//' @param x1         (n1 x p) numerator sample
//' @param x2         (n2 x p) denominator sample
//' @param centers    (M x p)  kernel centres
//' @param sigma_init (p)      starting bandwidths
//' @param n_y_cols   int      number of Y columns (get sigma_Y)
//' @param fold       int      CV folds
//' @param n_cd_iter  int      coordinate descent iterations
//' @return List with kernel_weights (vec), sigma (vec), centers (mat)
// [[Rcpp::export]]
Rcpp::List kliep_bw_cpp(const arma::mat& x1,
                         const arma::mat& x2,
                         const arma::mat& centers,
                         arma::vec sigma_init,
                         int n_y_cols,
                         int fold = 5,
                         int n_cd_iter = 3)
{
  const int p = (int)x1.n_cols;
  arma::vec sigma_vec = sigma_init;

  // build CV splits (1-indexed, 1..fold)
  int nx1 = (int)x1.n_rows;
  arma::ivec cv_split(nx1);
  arma::uvec perm = arma::randperm(nx1);
  for (int i = 0; i < nx1; i++) {
    cv_split[perm[i]] = (i % fold) + 1;
  }

  // coordinate indices for Y and X
  arma::uvec idx_Y(n_y_cols);
  for (int j = 0; j < n_y_cols; j++) idx_Y[j] = j;
  arma::uvec idx_X(p - n_y_cols);
  for (int j = n_y_cols; j < p; j++) idx_X[j - n_y_cols] = j;

  // coordinate descent
  double best_score = -1e30;
  for (int cd = 0; cd < n_cd_iter; cd++) {
    // optimise sigma_Y
    search_one_coord(sigma_vec[0], x1, x2, centers, sigma_vec,
                     idx_Y, cv_split, fold);

    // optimise sigma_X (if X columns exist)
    if (p > n_y_cols) {
      search_one_coord(sigma_vec[n_y_cols], x1, x2, centers, sigma_vec,
                       idx_X, cv_split, fold);
    }

    // check convergence
    arma::vec inv_2s2(p);
    for (int j = 0; j < p; j++) {
      inv_2s2[j] = 1.0 / (2.0 * sigma_vec[j] * sigma_vec[j]);
    }
    double cur = kliep_cv_score(x1, x2, centers, inv_2s2, cv_split, fold);
    if (cur <= best_score) break;
    best_score = cur;
  }

  // final weight optimisation with chosen sigma
  arma::vec inv_2s2(p);
  for (int j = 0; j < p; j++) {
    inv_2s2[j] = 1.0 / (2.0 * sigma_vec[j] * sigma_vec[j]);
  }
  arma::mat phi_x1 = compute_phi(x1, centers, inv_2s2);
  arma::mat phi_x2 = compute_phi(x2, centers, inv_2s2);
  arma::vec mean_phi_x2 = arma::mean(phi_x2, 0).t();
  arma::vec alpha = kliep_optimize(phi_x1, mean_phi_x2);

  return Rcpp::List::create(
    Rcpp::Named("kernel_weights") = Rcpp::wrap(alpha),
    Rcpp::Named("sigma")          = Rcpp::wrap(sigma_vec),
    Rcpp::Named("centers")        = Rcpp::wrap(centers)
  );
}


// ================================================================
//  EIF grid utilities
// ================================================================

//' Combined rowSums for E_alpha and E_Yalpha in one pass
//'
//' For each row i computes
//'   E_alpha[i]  = sum_j OR[i,j] * Cond[i,j]
//'   E_Yalpha[i] = sum_j Y_grid[j] * OR[i,j] * Cond[i,j]
//' in a single C++ loop, avoiding the creation of temporary N x M matrices.
//'
//' @param OR     (N x M) odds-ratio grid matrix
//' @param Cond   (N x M) conditional density grid  (already scaled by dY)
//' @param Y_grid (M)     outcome grid values
//' @return List with E_alpha (N) and E_Yalpha (N)
// [[Rcpp::export]]
Rcpp::List rowsums_grid_cpp(const arma::mat& OR,
                             const arma::mat& Cond,
                             const arma::vec& Y_grid)
{
  const int N = (int)OR.n_rows;
  const int M = (int)OR.n_cols;
  arma::vec E_alpha(N), E_Yalpha(N);

#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < N; i++) {
    double ea = 0.0, eya = 0.0;
    for (int j = 0; j < M; j++) {
      double v = OR(i, j) * Cond(i, j);
      ea  += v;
      eya += Y_grid[j] * v;
    }
    E_alpha[i]  = ea;
    E_Yalpha[i] = eya;
  }
  return Rcpp::List::create(Rcpp::Named("E_alpha")  = E_alpha,
                             Rcpp::Named("E_Yalpha") = E_Yalpha);
}


//' Sensitivity-scaled rowSums in one pass (no temporary N x M matrices)
//'
//' For each row i and column j applies the sensitivity scale factor
//'   scale[i,j] = exp(+half_log_Gamma)  if Y_grid[j] > mu[i]  (UB direction)
//'              = exp(-half_log_Gamma)   otherwise
//' directly during accumulation, avoiding construction of the scale and
//' scaled-OR matrices (saves 3-5 temporary N x M allocations per call).
//'
//' @param OR_grid       (N x M) base odds-ratio grid (alpha_0)
//' @param Cond          (N x M) conditional density grid
//' @param Y_grid        (M)     outcome grid values
//' @param mu_vec        (N)     per-observation baseline counterfactual mean
//' @param half_log_Gamma  scalar = log(Gamma) / 2
//' @param is_UB         bool    TRUE for upper-bound, FALSE for lower-bound
//' @return List with E_alpha (N) and E_Yalpha (N)
// [[Rcpp::export]]
Rcpp::List sens_rowsums_cpp(const arma::mat& OR_grid,
                             const arma::mat& Cond,
                             const arma::vec& Y_grid,
                             const arma::vec& mu_vec,
                             double           half_log_Gamma,
                             bool             is_UB)
{
  const int    N      = (int)OR_grid.n_rows;
  const int    M      = (int)OR_grid.n_cols;
  const double g_up   = std::exp( half_log_Gamma);
  const double g_down = std::exp(-half_log_Gamma);
  arma::vec E_alpha(N), E_Yalpha(N);

  // Branch outside the parallel loop so the condition is not re-evaluated N times
  if (is_UB) {
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < N; i++) {
      const double mu_i = mu_vec[i];
      double ea = 0.0, eya = 0.0;
      for (int j = 0; j < M; j++) {
        double scale = (Y_grid[j] > mu_i) ? g_up : g_down;
        double v     = OR_grid(i, j) * Cond(i, j) * scale;
        ea  += v;
        eya += Y_grid[j] * v;
      }
      E_alpha[i]  = ea;
      E_Yalpha[i] = eya;
    }
  } else {
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < N; i++) {
      const double mu_i = mu_vec[i];
      double ea = 0.0, eya = 0.0;
      for (int j = 0; j < M; j++) {
        double scale = (Y_grid[j] > mu_i) ? g_down : g_up;
        double v     = OR_grid(i, j) * Cond(i, j) * scale;
        ea  += v;
        eya += Y_grid[j] * v;
      }
      E_alpha[i]  = ea;
      E_Yalpha[i] = eya;
    }
  }
  return Rcpp::List::create(Rcpp::Named("E_alpha")  = E_alpha,
                             Rcpp::Named("E_Yalpha") = E_Yalpha);
}


//' Multiplier bootstrap SD (memory-efficient, uses R RNG)
//'
//' Computes sd(colMeans(matrix(rnorm(N*NumBoot)+1, N, NumBoot) * V))
//' using a single BLAS matrix-vector product instead of materialising
//' the full N x NumBoot matrix in R.  Uses R's RNG so set.seed() in R
//' controls reproducibility.
//'
//' Math:  colMeans(BootMat * V)[b]
//'          = mean_V + (1/N) * sum_i V[i] * Z[i,b]   (Z ~ N(0,1))
//'          = mean_V + (V^T Z)_b / N                  (BLAS DGEMV)
//'
//' @param V       length-N EIF vector
//' @param NumBoot number of bootstrap replicates
//' @return scalar standard deviation (same as sd(Mboot(V, N, NumBoot)))
// [[Rcpp::export]]
double mboot_sd_cpp(const arma::vec& V, int NumBoot)
{
  const int    N      = (int)V.n_elem;
  const double mean_V = arma::mean(V);

  // Draw N*NumBoot N(0,1) via R's RNG — respects set.seed() from R
  Rcpp::NumericVector z_r = Rcpp::rnorm(N * NumBoot, 0.0, 1.0);

  // Zero-copy view as column-major N x NumBoot Armadillo matrix
  const arma::mat Z(z_r.begin(), N, NumBoot, /*copy=*/false);

  // col_means[b] = mean_V + (V^T Z_col_b) / N  — single BLAS call
  const arma::rowvec col_means = V.t() * Z / N + mean_V;

  return arma::stddev(col_means);   // sd with 1/(N-1) denominator, same as R's sd()
}