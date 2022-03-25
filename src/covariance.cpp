#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(Rcpp)]]

// [[Rcpp::export]]
double covariance(NumericVector x, NumericVector y){
  double covar = 0.0;
  double mean_x = 0.0;
  double mean_y = 0.0;
  int n = x.length() - 1;
  
  mean_x = Rcpp::mean(x);
  mean_y = Rcpp::mean(y);
  
  for(int i = 0; i <= n; i++) {
    covar = covar + ((x(i) - mean_x) * (y(i) - mean_y));
  }
  covar = covar / n;
  
  return covar;
}