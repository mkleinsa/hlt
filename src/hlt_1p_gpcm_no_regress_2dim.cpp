#define RCPPDIST_DONT_USE_ARMA
#include <RcppDist.h>
using namespace Rcpp;
// [[Rcpp::depends(Rcpp, RcppDist)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include "abs2.h"

// [[Rcpp::export]]
double lgp1PNR2D(IntegerMatrix & x,
                 NumericVector lambda,
                 int nT,
                 int n,
                 int J,
                 NumericVector tJ,
                 int nDmax,
                 NumericVector lJ,
                 NumericVector theta,
                 NumericVector d,
                 double eps) {
  
  double llk = 0;
  d.attr("dim") = Dimension(nDmax, J);
  NumericMatrix dne = as<NumericMatrix>(d);
  
  theta.attr("dim") = Dimension(n, nT);
  NumericMatrix theta_mat = as<NumericMatrix>(theta);
  NumericMatrix d2(nDmax + 1, J);
  for(int i = 1; i <= nDmax; i++) {
    for(int j = 0; j < J; j++) {
      if(i <= lJ(j)) {
        d2(i, j) = dne(i - 1, j);
        double dt = d2(i, j);
        llk = llk + std::log(R::dnorm(dt, 0.0, 10.0, false) + eps);
      }
    }
  }
  
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < J; j++) {
      NumericVector r = cumsum((theta_mat(i, tJ(j)) - d2(_, j)));
      NumericVector er = exp(r);
      NumericVector pr = er / Rcpp::sum(er);
      int xval = x(i, j);
      double lpr = std::log(pr(xval) + eps);
      llk += lpr;
    }
  }
  
  for(int i = 0; i < n; i++) {
    double thetagt = theta_mat(i, nT - 1);
    for(int tt = 0; tt < (nT - 1); tt++) {
      double thetait = theta_mat(i, tt);
      double lambdat = lambda(tt);
      llk = llk + std::log(R::dnorm(thetait - (thetagt * lambdat), 0.0, 1.0, false) + eps);
    }
  }
  
  for(int i = 0; i < n; i++) {
    double thetagt = theta_mat(i, nT - 1);
    llk = llk + std::log(R::dnorm(thetagt, 0.0, 1.0, false) + eps);
  }
  
  for(int l = 0; l < nT - 1; l++) {
    double lambdal = lambda(l);
    llk = llk + std::log(R::dnorm(lambdal, 0.0, 1.0, false) + eps);
  }
  
  return llk;
}

// [[Rcpp::export]]
double lt1PNR2D(IntegerMatrix & x,
                int iter,
                int burn,
                double delta,
                NumericMatrix & post,
                NumericVector ix,
                NumericVector ixe,
                int npar,
                int n,
                int J,
                int nDmax,
                NumericVector lJ,
                int nT,
                NumericVector tJ,
                NumericVector & accept,
                double eps,
                bool display_progress = true) {
  
  Progress p(iter, display_progress);
  
  NumericVector oldpars = post(0, _ );
  oldpars[Range(ix(0), ix(0))] = oldpars[Range(ix(0) - 1, ix(0) - 1)];
  
  for(int it = 1; it < iter; it++) {
    NumericVector prop = Rcpp::rnorm(npar, 0.0, delta);
    NumericVector newpars = oldpars + prop;
    
    newpars[Range(ix(0), ix(0))] = newpars[Range(ix(0) - 1, ix(0) - 1)];
    
    double numer = lgp1PNR2D(x,
                             newpars[Range(ix(0) - 1, ixe(0) - 1)],
                             nT,
                             n,
                             J,
                             tJ,
                             nDmax,
                             lJ,
                             newpars[Range(ix(1) - 1, ixe(1) - 1)],
                             newpars[Range(ix(2) - 1, ixe(2) - 1)],
                             eps);
    
    double denom = lgp1PNR2D(x,
                             oldpars[Range(ix(0) - 1, ixe(0) - 1)],
                             nT,
                             n,
                             J,
                             tJ,
                             nDmax,
                             lJ,
                             oldpars[Range(ix(1) - 1, ixe(1) - 1)],
                             oldpars[Range(ix(2) - 1, ixe(2) - 1)],
                             eps);
    
    double acceptp = std::exp(numer - denom);
    double acceptit = (acceptp > R::runif(0.0, 1.0));
    
    if(acceptit == true) {
      oldpars = newpars;
      if(it >= burn) {
        post(it - burn, _ ) = newpars;
      }
      accept[it] = 1;
    } else {
      if(it >= burn) {
        post(it - burn, _ ) = oldpars;
      }
      accept[it] = 0;
    }
    
    p.increment();
  }
  
  return 1.0;
}
