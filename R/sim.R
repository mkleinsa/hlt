#' Simulate the HLT model
#' 
#' @param n number of observations
#' @param ntheta number first-level of latent dimensions
#' @param lambda latent factor coefficients
#' @param tJ number of questions
#' @param dL number of levels for each question
#' @param mua mean of the alphas
#' @param mud mean of the deltas
#' @param siga standard deviation of the alphas
#' @param sigd standard deviation of the deltas
#' @param regression TRUE/FALSE. Simulate regression model?
#' @param z regression design matrix.
#' @param nB number of regression parameters. nB = ncol(z).
#' @param beta what value to set the regression parameters.
#' 
#' @importFrom truncnorm rtruncnorm
#' @export
hltsim = function(n, ntheta, lambda, tJ, dL, mua, mud, siga, sigd, 
                  regression = FALSE, z = NULL, nB, beta = NULL) {
  nT = ntheta + 1
  theta = matrix(0, n, nT)
  theta[, nT] = seq(-3, 3, length.out = n)
  s.beta = beta #runif(nB, -2, 2)
  if(regression == TRUE) {
    theta[, nT] = theta[, nT] + z %*% s.beta + rnorm(n, 0, 0.1)
  }
  for(i in 1:(nT - 1)) {
    theta[, i] = lambda[i] * theta[, nT] + rnorm(n, 0, 1) #sqrt(1 - (lambda[i] ^ 2))
  }
  s.lam.cor = cor.theta(theta)
  J = length(tJ)
  nD = dL - 1
  s.theta = theta
  s.mua = mua
  s.mud = mud
  s.siga = siga
  s.sigd = sigd
  if(!is.null(s.mua)) {
    s.alpha = rtruncnorm(J, a = 0, mean = s.mua, sd = s.siga) #rnorm(J, s.mua, s.sigsqa)
  } else {
    s.alpha = 0
  }
  s.delta = mapply(1:J, FUN = function(x) {sort(rnorm(dL, s.mud, s.sigd))}, SIMPLIFY = "matrix")
  s.delta[1, ] = rep(0, J)
  s.lambda = lambda
  x = matrix(0, n, J)
  for (i in 1:n) {
    for (j in 1:J) {
      if(!is.null(s.mua)) {
        exp_part = exp(cumsum(s.alpha[j] * (s.theta[i, tJ[j] + 1] - s.delta[, j])))
      } else {
        exp_part = exp(cumsum((s.theta[i, tJ[j] + 1] - s.delta[, j])))
      }
      x[i, j] = sample(1:(dL) - 1, size = 1, prob = exp_part / sum(exp_part))
    }
  }
  return(list(x = x, 
              z = z,
              theta = theta,
              namesx = paste0("x", 1:J),
              s.beta = s.beta,
              s.lam.cor = s.lam.cor,
              s.mua = s.mua,
              s.mud = s.mud,
              s.siga = s.siga,
              s.sigd = s.sigd,
              s.alpha = s.alpha,
              s.delta = s.delta,
              s.lambda = s.lambda))
}
