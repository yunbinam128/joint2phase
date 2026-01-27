#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix update_p_vl_cpp(NumericMatrix prob_y0_v, NumericMatrix B0, 
                              NumericMatrix p_vl_curr, NumericMatrix p_vl_R1, 
                              int max_iter, double tol) {
  int n0 = prob_y0_v.nrow();
  int d_size = prob_y0_v.ncol();
  int s_n = B0.ncol();
  
  NumericMatrix p_vl = clone(p_vl_curr);
  
  for (int iter = 0; iter < max_iter; iter++) {
    NumericMatrix p_vl_old = clone(p_vl);
    NumericMatrix p_vl_num_s0(d_size, s_n);
    
    // Combined E and M steps to avoid 3D array overhead
    for (int i = 0; i < n0; i++) {
      double row_denom = 0.0;
      // Compute denominator for observation i
      for (int l = 0; l < s_n; l++) {
        for (int v = 0; v < d_size; v++) {
          row_denom += prob_y0_v(i, v) * B0(i, l) * p_vl(v, l);
        }
      }
      row_denom += 1e-12;
      
      // Update numerator sums
      for (int l = 0; l < s_n; l++) {
        for (int v = 0; v < d_size; v++) {
          p_vl_num_s0(v, l) += (prob_y0_v(i, v) * B0(i, l) * p_vl(v, l)) / row_denom;
        }
      }
    }
    
    // Update and normalize p_vl
    for (int l = 0; l < s_n; l++) {
      double col_sum = 0.0;
      for (int v = 0; v < d_size; v++) {
        p_vl(v, l) = p_vl_R1(v, l) + p_vl_num_s0(v, l);
        col_sum += p_vl(v, l);
      }
      if (col_sum == 0) col_sum = 1.0;
      for (int v = 0; v < d_size; v++) {
        p_vl(v, l) /= col_sum;
      }
    }
    
    // Check convergence
    double max_diff = 0;
    for (int j = 0; j < p_vl.length(); j++) {
      max_diff = std::max(max_diff, std::abs(p_vl[j] - p_vl_old[j]));
    }
    if (max_diff < tol) break;
  }
  return p_vl;
}