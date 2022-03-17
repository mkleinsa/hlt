#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(Rcpp)]]

// [[Rcpp::export]]
double abs2(double x) {
  double nx = 0.001;
  if(x < 0.0) {
    return nx;
  } else {
    return x;
  }
}