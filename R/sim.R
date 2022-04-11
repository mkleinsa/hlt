#' Simulate the HLT model
#' 
#' @param n number of observations
#' @param ntheta number first-level of latent dimensions
#' @param lambda latent factor coefficients
#' @param id number of questions
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
#' 
#' @examples 
#' 
#' # PCM, No regression, 4 dimensions
#' xdat = hltsim(type = "1p", n = 250, ntheta = 4, lambda = c(0.7, 0.2, 1.2, 0.4), 
#'               id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#'                      3,3,3,3,3,3), dL = 2)
#' apply(xdat$x, 2, table)
#' mod = hlt(xdat$x, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#'           3,3,3,3,3,3), iter = 200, burn = 100, delta = 0.013)
#' 
hltsim = function(type, n, ntheta, lambda, id, dL, nB, beta = NULL) {
  
  if(!is.null(beta)) {
    regression = TRUE
  } else {
    regression = FALSE
  }
  
  nT = ntheta + 1
  
  theta = matrix(0, n, nT)
  
  theta[, nT] = seq(-3, 3, length.out = n)
  
  s.beta = beta
  
  if(regression == TRUE) {
    z = matrix(sample(c(0, 1), size = n * nB, replace = TRUE), nrow = n, 
               ncol = nB)
    theta[, nT] = theta[, nT] + z %*% s.beta + rnorm(n, 0, 0.1)
  }
  
  for(i in 1:(nT - 1)) {
    theta[, i] = lambda[i] * theta[, nT] + rnorm(n, 0, 1)
  }
  
  s.lam.cor = cor.theta(theta)
  
  J = length(id)
  
  nD = dL - 1
  
  s.theta = theta
  
  s.alpha = rtruncnorm(J, a = 0, mean = 1, sd = 0.2)
  
  s.delta = mapply(1:J, FUN = function(x) {sort(rnorm(dL, 1, 0.2))}, 
                   SIMPLIFY = "matrix")
  
  s.delta[1, ] = rep(0, J)
  
  s.lambda = lambda
  
  x = matrix(0, n, J)
  for (i in 1:n) {
    for (j in 1:J) {
      if(type == "2p") {
        exp_part = exp(cumsum((s.alpha[j] * (s.theta[i, id[j] + 1])) - s.delta[, j]))
      } else if(type == "1p") {
        exp_part = exp(cumsum((s.theta[i, id[j] + 1]) - s.delta[, j]))
      }
      x[i, j] = sample(1:(dL) - 1, size = 1, prob = exp_part / sum(exp_part))
    }
  }
  
  if(regression == TRUE) {
    return(list(x = x, 
                z = z,
                id = id,
                theta = theta,
                namesx = paste0("x", 1:J),
                s.beta = s.beta,
                s.lam.cor = s.lam.cor,
                s.alpha = s.alpha,
                s.delta = s.delta,
                s.lambda = s.lambda))
  } else {
    return(list(x = x, 
                theta = theta,
                id = id,
                namesx = paste0("x", 1:J),
                s.lam.cor = s.lam.cor,
                s.alpha = s.alpha,
                s.delta = s.delta,
                s.lambda = s.lambda))
  }
}
