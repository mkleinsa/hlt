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
#' 
#' @examples 
#' xdat = hltsim(n = 100, ntheta = 3, lambda = c(0.8, 0.4, 0.3), 
#'               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2),
#'               dL = 5, mua = 1, mud = 1.4, siga = 0.8, sigd = 1,
#'               regression = FALSE, z = NULL)
#' apply(xdat$x, 2, table)
#'               
#' nB = 1
#' n = 100
#' z = matrix(sample(0:1, n, replace = TRUE), nrow = n, ncol = nB)
#' xdat = hltsim(n = 100, ntheta = 3, lambda = c(0.8, 0.4, 0.3), 
#'               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2),
#'               dL = 5, mua = 1, mud = 1.4, siga = 0.8, sigd = 1,
#'               regression = TRUE, z = z, nB = nB)
#' apply(xdat$x, 2, table)
#' lm(xdat$theta[,4] ~ xdat$z)
#' xdat$s.beta
#' 
#' @importFrom truncnorm rtruncnorm
hltsim = function(n, ntheta, lambda, tJ, dL, mua, mud, siga, sigd, 
                  regression = FALSE, z = NULL, nB, beta) {
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
  s.alpha = rtruncnorm(J, a = 0, mean = s.mua, sd = s.siga) #rnorm(J, s.mua, s.sigsqa)
  s.delta = mapply(1:J, FUN = function(x) {sort(rnorm(dL, s.mud, s.sigd))}, SIMPLIFY = "matrix")
  s.delta[1, ] = rep(0, J)
  s.lambda = lambda
  x = matrix(0, n, J)
  for (i in 1:n) {
    for (j in 1:J) {
      exp_part = exp(cumsum(s.alpha[j] * (s.theta[i, tJ[j] + 1] - s.delta[, j])))
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
