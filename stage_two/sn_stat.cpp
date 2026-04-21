#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// 1. 计算单变量全局 SN 统计量 (用于展示 SN+BS 的失效)
// 严格对应论文公式 (3) 的 Vn(k) 结构
// [[Rcpp::export]]
NumericVector calc_global_sn_univariate(arma::vec Y) {
    int n = Y.n_elem;
    NumericVector Tn(n);
    arma::vec S = arma::cumsum(Y);
    double total_sum = S(n-1);

    int min_size = 30; // 边缘截断
    for (int k = min_size; k < n - min_size; k++) { 
        double k_val = (double)k + 1.0;
        
        // Dn(k) - 公式 (3) 分子
        double mean_left = S(k) / k_val;
        double mean_right = (total_sum - S(k)) / (n - k_val);
        double Dn_k = (k_val * (n - k_val) / pow(n, 1.5)) * (mean_left - mean_right);
        
        // Ln(k) - 公式 (3) 左侧方差
        double L_sum = 0;
        for (int i = 1; i <= (int)k_val; i++) {
            double i_val = (double)i;
            double term = (S(i-1)/i_val - mean_left);
            L_sum += (pow(i_val, 2) * pow(k_val - i_val, 2)) * pow(term, 2);
        }
        double Ln = L_sum / (pow(n, 2) * pow(k_val, 2));
        
        // Rn(k) - 公式 (3) 右侧方差
        double R_sum = 0;
        for (int i = (int)k_val + 1; i <= n; i++) {
            double i_val = (double)i;
            double term = ((total_sum - S(i-2))/(n - i_val + 1.0) - mean_right);
            R_sum += (pow(n - i_val + 1.0, 2) * pow(i_val - k_val - 1.0, 2)) * pow(term, 2);
        }
        double Rn = R_sum / (pow(n, 2) * pow(n - k_val, 2));
        
        double Vn = Ln + Rn;
        if (Vn > 1e-10) Tn[k] = pow(Dn_k, 2) / Vn;
    }
    return Tn;
}

// 2. 多变量 CUSUM 统计量 (用于加速 CUSUM_Multi 和 AHHR)
// [[Rcpp::export]]
NumericVector calc_cusum_multi(arma::mat Y) {
    int n = Y.n_rows;
    NumericVector stats(n);
    arma::mat S = arma::cumsum(Y, 0); 
    
    int min_size = 30;
    for(int k = min_size; k < n - min_size; k++) {
        double k_val = (double)k + 1.0;
        arma::rowvec mean_left = S.row(k) / k_val;
        arma::rowvec mean_right = (S.row(n-1) - S.row(k)) / (n - k_val);
        arma::rowvec diff = mean_left - mean_right;
        
        // 计算 L2 范数作为多元跳变度量
        stats[k] = std::sqrt(k_val * (n - k_val) / n) * arma::norm(diff, 2);
    }
    return stats;
}