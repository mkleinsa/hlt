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
           int nD,
           NumericVector theta,
           NumericVector d,
           double mud,
           double sigd,
           NumericVector a,
           NumericVector beta,
           double eps,
           double mud_prior_mean = 0.0,
           double mud_prior_stdev = 10.0,
           double sigd_prior_max = 2.0,
           double lambdal_prior_min = 0.0,
           double lambdal_prior_max = 2.0) {

  double llk = 0;

  d.attr("dim") = Dimension(nD, J);
  NumericMatrix dne = as<NumericMatrix>(d);

  theta.attr("dim") = Dimension(n, nT);
  NumericMatrix theta_mat = as<NumericMatrix>(theta);

  NumericMatrix d2(nD + 1, J);
  for(int i = 1; i <= nD; i++) {
    for(int j = 0; j < J; j++) {
      d2(i, j) = dne(i - 1, j);
      double dt = d2(i, j);
      llk = llk + std::log(R::dnorm((dt - mud) / sigd, 0.0, 1.0, false) + eps);
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
    // llk = llk + std::log(R::dunif(lambdal, lambdal_prior_min, lambdal_prior_max,
    //                               false) + eps);
    llk = llk + std::log(d_truncnorm(lambdal, 0.0, 2.0, 0.0, 1000.0) + eps);
  }

  for(int j = 0; j < J; j++) {
    double at = a(j);
    llk = llk + std::log(d_truncnorm(at, 0.0, 10.0, 0.0, 1000.0) + eps);
  }
  
  for(int zz = 0; zz < nB; zz++) {
    double betaz = beta(zz);
    llk = llk + std::log(R::dnorm(betaz, 0.0, 10.0, false) + eps);
  }
  
  llk = llk + std::log(R::dnorm(mud, mud_prior_mean, mud_prior_stdev, false) + eps);
  llk = llk + std::log(R::dunif(sigd, 0.0, sigd_prior_max, false) + eps);

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
          int nD,
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

    for(int q = ix(0) - 1; q < ixe(0); q ++){
      oldpars(q) = abs2(oldpars(q));
      newpars(q) = abs2(newpars(q));
    }
    for(int q = ix(4) - 1; q < ixe(4); q ++){
      oldpars(q) = abs2(oldpars(q));
      newpars(q) = abs2(newpars(q));
    }
    for(int q = ix(5) - 1; q < ixe(5); q ++){
      oldpars(q) = abs2(oldpars(q));
      newpars(q) = abs2(newpars(q));
    }

    double numer = lgp(x,
                       z,
                       newpars[Range(ix(0), ixe(0))],
                       nB,
                       nT,
                       n,
                       J,
                       tJ,
                       nD,
                       newpars[Range(ix(1), ixe(1))],
                       newpars[Range(ix(2), ixe(2))],
                       newpars[Range(ix(3), ixe(3))][0],
                       newpars[Range(ix(4), ixe(4))][0],
                       newpars[Range(ix(5), ixe(5))],
                       newpars[Range(ix(6), ixe(6))],
                       eps);

    double denom = lgp(x,
                       z,
                       oldpars[Range(ix(0), ixe(0))],
                       nB,
                       nT,
                       n,
                       J,
                       tJ,
                       nD,
                       oldpars[Range(ix(1), ixe(1))],
                       oldpars[Range(ix(2), ixe(2))],
                       oldpars[Range(ix(3), ixe(3))][0],
                       oldpars[Range(ix(4), ixe(4))][0],
                       oldpars[Range(ix(5), ixe(5))],
                       newpars[Range(ix(6), ixe(6))],
                       eps);

    double acceptp = std::exp(numer - denom);
    double acceptit = (acceptp > R::runif(0.0, 1.0));
    
    if(acceptit == true) {
      oldpars = newpars;
      if(it > burn) {
        post(it - burn, _ ) = newpars;
      }
      accept[it] = 1;
    } else {
      if(it > burn) {
        post(it - burn, _ ) = oldpars;
      }
      accept[it] = 0;
    }
    
    p.increment();
  }
  
  return 1.0;
}
