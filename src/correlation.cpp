#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(Rcpp)]]

// [[Rcpp::export]]
double correlation(NumericVector x, NumericVector y){
  double covar = 0.0;
  double var_x = 0.0;
  double var_y = 0.0;
  double corr = 0.0;
  double mean_x = 0.0;
  double mean_y = 0.0;
  int n = x.length() - 1;
  
  mean_x = Rcpp::mean(x);
  mean_y = Rcpp::mean(y);
  
  for(int i = 0; i <= n; i++) {
    covar = covar + ((x(i) - mean_x) * (y(i) - mean_y));
    var_x = var_x + ((x(i) - mean_x) * (x(i) - mean_x));
    var_y = var_y + ((y(i) - mean_y) * (y(i) - mean_y));
  }
  covar = covar / n;
  var_x = var_x / n;
  var_y = var_y / n;
  
  corr = covar / (std::sqrt(var_x) * std::sqrt(var_y));
  
  return corr;
}