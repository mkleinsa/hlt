#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(Rcpp)]]
#include "correlation.h"

// [[Rcpp::export]]
void standardize_lambda(Rcpp::NumericMatrix & post, int start_idx, int end_idx,
                        int nT, int n, Rcpp::NumericMatrix & corr_theta) {
  
  NumericMatrix theta_draws = post( _ , Rcpp::Range(start_idx, end_idx));
  int nr = post.nrow();
  
  for(int draw = 0; draw < nr; draw++) {
    NumericVector theta = theta_draws(draw, _ );
    theta.attr("dim") = Dimension(n, nT);
    NumericMatrix theta_mat = as<NumericMatrix>(theta);
    NumericVector theta_g = theta_mat( _ , nT - 1);
    for(int cr = 0; cr < (nT - 1); cr++) {
      corr_theta(draw, cr) = correlation(theta_mat( _ , cr), theta_g);
    }
  }
  
}