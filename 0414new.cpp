#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// 计算多维特征差的 L2 范数平方 (核心检验统计量 S)
double calc_stat(const mat& Z, int L) {
  rowvec mean_before = mean(Z.rows(0, L - 1), 0);
  rowvec mean_after  = mean(Z.rows(L, 2 * L - 1), 0);
  return accu(square(mean_after - mean_before));
}

// [[Rcpp::export]]
double boot_test_cpp(mat Z, int L, int block_size, int M) {
  int n = 2 * L;
  int d = Z.n_cols;
  
  mat Z_norm = Z;
  for(int j = 0; j < d; ++j) {
    double col_sd = stddev(Z.col(j));
    if(col_sd > 1e-6) {
      Z_norm.col(j) /= col_sd;
    }
  }
  
  // 1. 计算观测到的真实统计量
  double S_obs = calc_stat(Z_norm, L);
  
  // 2. 强加原假设 (Null-adjustment): 消除 t > L 的均值差
  rowvec mean_before = mean(Z_norm.rows(0, L - 1), 0);
  rowvec mean_after  = mean(Z_norm.rows(L, n - 1), 0);
  rowvec diff = mean_after - mean_before;
  
  mat Z_tilde = Z_norm;
  for(int t = L; t < n; ++t) {
    Z_tilde.row(t) -= diff;
  }
  
  // 3. Moving Block Bootstrap (MBB) 抽样与检验
  int n_blocks = std::ceil((double)n / block_size);
  int exceed_count = 0;
  
  for(int m = 0; m < M; ++m) {
    mat Z_star(n, d, fill::zeros);
    int current_idx = 0;
    
    for(int b = 0; b < n_blocks; ++b) {
      int start_idx = R::runif(0, n - block_size + 1); 
      for(int k = 0; k < block_size; ++k) {
        if (current_idx < n) {
          Z_star.row(current_idx) = Z_tilde.row(start_idx + k);
          current_idx++;
        }
      }
    }
    
    double S_star = calc_stat(Z_star, L);
    if(S_star > S_obs) {
      exceed_count++;
    }
  }
  
  return (double)exceed_count / M;
}