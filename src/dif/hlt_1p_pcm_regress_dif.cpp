#define RCPPDIST_DONT_USE_ARMA
#include <RcppDist.h>
using namespace Rcpp;
// [[Rcpp::depends(Rcpp, RcppDist)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::export]]
double lgp(IntegerMatrix & x,
           NumericMatrix & z,
           NumericVector lambda,
           int nB,
           int nT,
           int n,
           int J,
           NumericVector tJ,
           int nDmax,
           NumericVector lJ,
           NumericVector theta,
           NumericVector d,
           double mud,
           double sigd,
           NumericVector a,
           NumericVector beta,
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
      NumericVector r = cumsum((theta_mat(i, tJ(j)) - d2(_, j)) * a(j));
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
    NumericVector zi = z(i, _ );
    double zzsum = 0;
    for(int zz = 0; zz < nB; zz++) {
      zzsum += zi(zz) * beta(zz);
    }
    llk = llk + std::log(R::dnorm(thetagt - zzsum, 0.0, 1.0, false) + eps);
  }

  for(int l = 0; l < nT - 1; l++) {
    double lambdal = lambda(l);
    llk = llk + std::log(R::dnorm(lambdal, 0.0, 1.0, false) + eps);
  }

  for(int j = 0; j < J; j++) {
    double at = a(j);
    llk = llk + std::log(d_truncnorm(at, 0.0, 10.0, 0.0, 1000.0) + eps);
  }
  
  for(int zz = 0; zz < nB; zz++) {
    double betaz = beta(zz);
    llk = llk + std::log(R::dnorm(betaz, 0.0, 10.0, false) + eps);
  }

  return llk;
}

// [[Rcpp::export]]
double abs2(double x) {
  double nx = 0.001;
  if(x < 0.0) {
    return nx;
  } else {
    return x;
  }
}

// [[Rcpp::export]]
double lt(IntegerMatrix & x,
          NumericMatrix & z,
          int iter,
          int burn,
          double delta,
          NumericMatrix & post,
          NumericVector ix,
          NumericVector ixe,
          int npar,
          int n,
          int nB,
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

  for(int it = 1; it < iter; it++) {
    NumericVector prop = Rcpp::rnorm(npar, 0.0, delta);
    NumericVector newpars = oldpars + prop;
    
    for(int q = ix(4) - 1; q < ixe(4); q ++){
      oldpars(q) = abs2(oldpars(q));
      newpars(q) = abs2(newpars(q));
    }
    for(int q = ix(5) - 1; q < ixe(5); q ++){
      oldpars(q) = abs2(oldpars(q));
      newpars(q) = abs2(newpars(q));
    }
    
    // NumericVector beta_old = oldpars[Range(ix(6) - 1, ixe(6) - 1)];
    // NumericVector beta_new = clone(beta_old);
    // int lb = beta_old.length();
    // int lthetag = (ixe(1) - ix(1)) / nT;
    // NumericVector thetag_old = oldpars[Range(ix(1) - 1 + lthetag - (lthetag / nT), ixe(1) - 1 + lthetag)];
    // 
    // double tau = 1.0 / 10.0;
    // for(int q = 0; q < lb; q++) {
    //   NumericVector zj = z( _ , q);
    //   NumericVector res = clone(thetag_old);
    //   NumericVector beta_nq = clone(beta_old);
    //   beta_nq[q] = 0.0;
    //   double sum_part = 0.0;
    //   double var_part = 0.0;
    //   for(int i = 0; i < lthetag; i++) {
    //     NumericVector zi = z( i , _ );
    //     double resid = sum(beta_nq * zi);
    //     double ziq = zi(q);
    //     sum_part = sum_part + (thetag_old(i) - resid) * ziq;
    //     var_part = var_part + (ziq * ziq);
    //   }
    //   beta_new(q) = R::rnorm(sum_part / (tau + var_part), std::sqrt(1.0 / (tau + var_part)));
    // }

    double numer = lgp(x,
                       z,
                       newpars[Range(ix(0) - 1, ixe(0) - 1)],
                       nB,
                       nT,
                       n,
                       J,
                       tJ,
                       nDmax,
                       lJ,
                       newpars[Range(ix(1) - 1, ixe(1) - 1)],
                       newpars[Range(ix(2) - 1, ixe(2) - 1)],
                       newpars[Range(ix(3) - 1, ixe(3) - 1)][0],
                       newpars[Range(ix(4) - 1, ixe(4) - 1)][0],
                       newpars[Range(ix(5) - 1, ixe(5) - 1)],
                       newpars[Range(ix(6) - 1, ixe(6) - 1)],
                       eps);

    double denom = lgp(x,
                       z,
                       oldpars[Range(ix(0) - 1, ixe(0) - 1)],
                       nB,
                       nT,
                       n,
                       J,
                       tJ,
                       nDmax,
                       lJ,
                       oldpars[Range(ix(1) - 1, ixe(1) - 1)],
                       oldpars[Range(ix(2) - 1, ixe(2) - 1)],
                       oldpars[Range(ix(3) - 1, ixe(3) - 1)][0],
                       oldpars[Range(ix(4) - 1, ixe(4) - 1)][0],
                       oldpars[Range(ix(5) - 1, ixe(5) - 1)],
                       oldpars[Range(ix(6) - 1, ixe(6) - 1)],
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
