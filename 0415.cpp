#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

static inline double calc_stat(const mat& Z, int L) {
  rowvec mean_before = mean(Z.rows(0, L - 1), 0);
  rowvec mean_after = mean(Z.rows(L, 2 * L - 1), 0);
  return accu(square(mean_after - mean_before));
}

static inline vec gen_dwb_weights(int n, double bandwidth) {
  double bw = std::max(1.0, bandwidth);
  double a = std::exp(-1.0 / bw);
  double s = std::sqrt(std::max(1e-12, 1.0 - a * a));

  vec xi(n, fill::zeros);
  xi(0) = R::rnorm(0.0, 1.0);
  for (int t = 1; t < n; ++t) {
    xi(t) = a * xi(t - 1) + s * R::rnorm(0.0, 1.0);
  }

  xi -= mean(xi);
  double sd_xi = stddev(xi);
  if (sd_xi > 1e-12) {
    xi /= sd_xi;
  }
  return xi;
}

// [[Rcpp::export]]
List dwb_test_cpp(mat Z, int L, int M = 500, double bandwidth = 10.0) {
  int n = Z.n_rows;
  int d = Z.n_cols;

  if (n != 2 * L) {
    stop("Z.n_rows must be exactly 2 * L.");
  }
  if (L < 2 || d < 1 || M < 10) {
    stop("Invalid inputs: require L >= 2, d >= 1, M >= 10.");
  }

  mat Z_norm = Z;
  for (int j = 0; j < d; ++j) {
    double col_sd = stddev(Z_norm.col(j));
    if (col_sd > 1e-8) {
      Z_norm.col(j) /= col_sd;
    }
  }

  double S_obs = calc_stat(Z_norm, L);

  rowvec mean_before = mean(Z_norm.rows(0, L - 1), 0);
  rowvec mean_after = mean(Z_norm.rows(L, n - 1), 0);
  rowvec diff = mean_after - mean_before;

  mat Z_tilde = Z_norm;
  for (int t = L; t < n; ++t) {
    Z_tilde.row(t) -= diff;
  }

  int exceed_count = 0;
  vec S_boot(M, fill::zeros);

  for (int m = 0; m < M; ++m) {
    vec xi = gen_dwb_weights(n, bandwidth);
    mat Z_star = Z_tilde;
    for (int t = 0; t < n; ++t) {
      Z_star.row(t) *= xi(t);
    }

    double S_star = calc_stat(Z_star, L);
    S_boot(m) = S_star;
    if (S_star >= S_obs) {
      exceed_count++;
    }
  }

  double p_value = static_cast<double>(exceed_count + 1) / static_cast<double>(M + 1);

  vec S_sorted = sort(S_boot);
  int idx95 = std::min(M - 1, std::max(0, static_cast<int>(std::floor(0.95 * M)) - 1));
  double S_q95 = S_sorted(idx95);
  double S_mean = mean(S_boot);

  return List::create(
    _["p_value"] = p_value,
    _["S_obs"] = S_obs,
    _["S_boot_q95"] = S_q95,
    _["S_boot_mean"] = S_mean
  );
}
