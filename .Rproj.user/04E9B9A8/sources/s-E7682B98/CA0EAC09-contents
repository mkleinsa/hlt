#' Higher-Order Item Response Theory (Latent Trait Theory) with Regression
#' 
#' Fit a higher-order item response theory model under the generalized 
#' partial credit measurement model. The goal is to explain multiple latent dimensions
#' by a single higher-order dimension. We extend this model with an option to perform regression 
#' on the general latent dimension.
#'
#' @param x matrix of item responses. Responses must be integers where the 
#' lowest value is 0 and the highest value is the maximum possible response
#' for the item with no gaps. If a question is asked with 5 possible responses,
#' then the possible values should be c(0,1,2,3,4). For binary items, use
#' c(0,1).
#' @param z centered numeric matrix of predictors for the latent regression. 
#' Default is `z = NULL` so that no regression is performed. All columns
#' of this matrix must be numeric. For binary items, use the values c(0,1).
#' For continuous items, center the values on the mean and divide by the 
#' standard deviation (normalized). For factors with more than two levels, 
#' recode into multiple columns of c(0,1).
#' @param id I.D. vector indexing first-order latent dimension membership
#' for each of the first-order latent dimensions. We index starting from zero,
#' not one. If there are three first-order .
#' latent dimensions with 5 questions per dimension, then the vector will look
#' like c(0,0,0,0,0,1,1,1,1,1,2,2,2,2,2).
#' @param iter number of total iterations.
#' @param burn number of burn in iterations.
#' @param delta tuning parameter for Metropolis-Hanstings algorithm. Alter 
#' delta until acceptance.ratio =~ 0.234.
#' @param type type of Partial Credit Model to fit. If the partial credit model
#' is desired (i.e. all alpha parameters = 1), then choose `type = "1p"`.
#' If the Generalized Parial Credit Model is desired, then choose
#' `type = "2p"`. The default is `type = "2p"`.
#' @param start starting values for the Metropolis-Hastings algorithm. 
#' Provide a `list` with the following named arguments:
#' `list(lambda = c(), theta = c(), delta = c(), alpha = c(), beta = c())`
#' \itemize{
#'   \item{lambda - }{vector of starting values for the latent factor loadings.}
#'   \item{theta - }{vector of starting values for the abilities.}
#'   \item{delta - }{vector of starting values for the difficulties.}
#'   \item{alpha - }{vector of starting values for the slope parameters.}
#'   \item{beta - }{vector of starting values for the latent regression parameters}
#' }
#' If you choose specify starting values, then the lengths of the starting value
#' vectors must match the number of parameters in the model.
#' 
#' @param progress boolean, show progress bar? Defaults to TRUE.
#'
#'
#' @return A matrix of posterior estimates. Rows are the draws and columns
#' are the named parameters.
#' 
#' @examples 
#' 
#' # Example 1: sumulated data
#' ntheta = 3
#' id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2)
#' lambda = c(0.8, 0.4, 0.3)
#' dL = 5
#' 
#' # small n test
#' set.seed(7794)
#' xdat = hltsim(n = 100, ntheta = ntheta, lambda = lambda, 
#'               tJ = id, dL = dL, mua = 1, mud = 1.4, siga = 0.8, sigd = 1)
#' x = xdat$x
#' mod = hlt(x, ntheta = 3, id = id, iter = 1000, delta = 0.03)
#' mod$accept.rate
#' post = mod$post
#' 
#' # increase n
#' set.seed(7794)
#' xdat = hltsim(n = 150, ntheta = ntheta, lambda = lambda, 
#'               tJ = id, dL = dL, mua = 1, mud = 1.4, siga = 0.8, sigd = 1,
#'               beta = 1)
#' x = xdat$x
#' mod = hlt(x, ntheta = 3, id = id, iter = 40000, burn = 30000, delta = 0.05)
#' mod$accept.rate
#' post = mod$post
#' 
#' # Regression example 
#' set.seed(2034)
#' nB = 1
#' n = 300
#' z = matrix(sample(0:1, n, replace = TRUE), nrow = n, ncol = nB)
#' xdat = hltsim(n = n, ntheta = 3, lambda = c(0.8, 0.4, 0.3), 
#'               tJ = id, dL = 4, mua = 1, mud = 1.4, siga = 0.8, sigd = 1,
#'               regression = TRUE, z = z, nB = nB, beta = 1)
#' apply(xdat$x, 2, table)
#' lm(xdat$theta[,4] ~ xdat$z)
#' xdat$s.beta
#' x = xdat$x
#' x[,1] = ifelse(x[,1] == 3, 2, x[,1])
#' x[,1] = ifelse(x[,1] == 2, 1, x[,1])
#' x[,4] = ifelse(x[,4] == 3, 2, x[,4])
#' x[,5] = ifelse(x[,5] == 3, 2, x[,5])
#' mod = hlt(x, z = z, id = id, iter = 1e5, delta = 0.017)
#' mod = hlt(x, z = z, id = id, iter = 1e6, burn = 9e5, delta = 0.017)
#' 
#' mod = hlt(x, z = z, id = id, iter = 2.2e6, burn = 2e6, delta = 0.013)
#' mod$accept.rate
#' post = mod$post
#' apply(post, 2, mean)
#' 
#' smy = function(x) {c(mean = mean(x), se = sd(x), quantile(x, probs = c(0.025, 0.5, 0.7, 0.975)))}
#' apply(post[, "beta1", drop = FALSE], 2, smy)
#' 
#' summary(mod, param = "all")
#' summary(mod, param = "beta")
#' summary(mod, param = "lambda")
#' summary(mod, param = "alpha")
#' summary(mod, param = "delta")
#' summary(mod, param = "theta", dimension = 1)
#' summary(mod, param = "theta", dimension = 2)
#' summary(mod, param = "theta", dimension = 3)
#' summary(mod, param = "theta", dimension = 4)
#' 
#' cor(as.vector(xdat$s.delta[-1,]), 
#'     as.vector(matrix(summary(mod, param = "delta")[1,], nrow = 3)))
#'     
#' cor(summary(mod, param = "alpha")[1,], xdat$s.alpha)
#' 
#' summary(mod, param = "cor.theta", cor.theta = c(1,2))
#' xdat$s.lam.cor
#' 
#' th = summary(mod, param = "theta", dimension = 1)
#' cor(th, xdat$theta[,1])
#' 
#' plot(mod, "lambda1", type = "trace")
#' plot(mod, "lambda2", type = "trace")
#' plot(mod, "lambda3", type = "trace")
#' plot(mod, "a1", type = "trace")
#' plot(mod, "d2", type = "trace")
#' plot(mod, "beta1", type = "trace")
#' 
#' plot(mod, 1, type = "icc")
#' plot(mod, 2, type = "icc")
#' plot(mod, 3, type = "icc")
#' plot(mod, 4, type = "icc")
#' plot(mod, 5, type = "icc")
#' plot(mod, 6, type = "icc")
#' plot(mod, 7, type = "icc", min = -10, max = 10)
#' 
#' # example 
#' data("asti")
#' x = as.matrix(asti[, 1:25]) - 1
#' z = asti[, 26:27]
#' z[, 1] = (z[, 1] == "students") * 1
#' z[, 2] = (z[, 2] == "male") * 1
#' z = as.matrix(z)
#' id = c(0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4)
#' mod = hlt(x, z = z, id = id, iter = 2.2e6, burn = 2e6, delta = 0.013)
#' 
#' 
#' ### Test each combination of input values ###
#' n = 200
#' nB = 5
#' z = matrix(sample(0:1, n, replace = TRUE), nrow = n, ncol = nB)
#' 
#' # PCM, No regression, 2 dimensions
#' xdat = hltsim(n = n, ntheta = 2, lambda = c(0.7, 0.7), 
#'               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
#'               dL = 2, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#'               regression = FALSE, z = NULL)
#' apply(xdat$x, 2, table)
#' mod = hlt(xdat$x, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1), 
#'           iter = 200, burn = 100, delta = 0.013)
#' 
#' # PCM, Regression, 2 dimensions
#' xdat = hltsim(n = n, ntheta = 2, lambda = c(0.7, 0.7), 
#'               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
#'               dL = 2, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#'               regression = FALSE, z = z)
#' apply(xdat$x, 2, table)
#' mod = hlt(xdat$x, z = z, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1), 
#'           iter = 200, burn = 100, delta = 0.013)
#' 
#' # PCM, No regression, 4 dimensions
#' xdat = hltsim(n = n, ntheta = 4, lambda = c(0.7, 0.2, 1.2, 0.4), 
#'               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#'                      3,3,3,3,3,3),
#'               dL = 2, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#'               regression = FALSE, z = NULL)
#' apply(xdat$x, 2, table)
#' mod = hlt(xdat$x, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#'           3,3,3,3,3,3), iter = 200, burn = 100, delta = 0.013)
#' 
#' # PCM, Regression, 4 dimensions
#' xdat = hltsim(n = n, ntheta = 4, lambda = c(0.7, 0.2, 1.2, 0.4), 
#'               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#'                      3,3,3,3,3,3),
#'               dL = 2, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#'               regression = FALSE, z = z)
#' apply(xdat$x, 2, table)
#' mod = hlt(xdat$x, z = z, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#'           3,3,3,3,3,3), iter = 200, burn = 100, delta = 0.013)
#' 
#' 
#' 
#' 
#' 
#' 
#' # GPCM, No regression, 2 dimensions
#' xdat = hltsim(n = n, ntheta = 2, lambda = c(0.7, 0.7), 
#'               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
#'               dL = 10, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#'               regression = FALSE, z = NULL)
#' apply(xdat$x, 2, table)
#' mod = hlt(xdat$x, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1), 
#'           iter = 200, burn = 100, delta = 0.013)
#' 
#' # GPCM, Regression, 2 dimensions
#' xdat = hltsim(n = n, ntheta = 2, lambda = c(0.7, 0.7), 
#'               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
#'               dL = 10, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#'               regression = FALSE, z = z)
#' apply(xdat$x, 2, table)
#' mod = hlt(xdat$x, z = z, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1), 
#'           iter = 200, burn = 100, delta = 0.013)
#' 
#' # GPCM, No regression, 4 dimensions
#' xdat = hltsim(n = n, ntheta = 4, lambda = c(0.7, 0.2, 1.2, 0.4), 
#'               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#'                      3,3,3,3,3,3),
#'               dL = 10, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#'               regression = FALSE, z = NULL)
#' apply(xdat$x, 2, table)
#' mod = hlt(xdat$x, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#'           3,3,3,3,3,3), iter = 200, burn = 100, delta = 0.013)
#' 
#' # GPCM, Regression, 4 dimensions
#' xdat = hltsim(n = n, ntheta = 4, lambda = c(0.7, 0.2, 1.2, 0.4), 
#'               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#'                      3,3,3,3,3,3),
#'               dL = 10, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#'               regression = FALSE, z = z)
#' apply(xdat$x, 2, table)
#' mod = hlt(xdat$x, z = z, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#'           3,3,3,3,3,3), iter = 200, burn = 100, delta = 0.013)
#' 
#' 
hlt = function(x, 
               z = NULL, 
               id, 
               iter, 
               burn = iter / 2, 
               delta,
               type = "2p",
               start = list(lambda = c(), theta = c(), delta = c(), 
                            alpha = c(), beta = c()),
               progress = TRUE) {
  
    if(!is.matrix(x)) {
      x = as.matrix(x)
    }
    
    ntheta = length(unique(id))
    
    if(ntheta < 2) {
      stop("Specified ntheta < 2. Must assume at least two latent dimensions
           to perform inference.")
    }
    
    if(ntheta == 2) {
      warning("Specified ntheta == 2. The lambda parameters for the two 
              specified dimensions are set to be equal, i.e. lambda1 == lambda2.
              This constraint is only required with < 3 dimensions.")
    }
    
    if(length(id) != ncol(x)) {
      stop("The id vector must match number of columns in x.")
    }
  
    n = nrow(x)
    nT = ntheta + 1
    J = ncol(x)
    lJ = apply(x, 2, function(x) {length(unique(x))})
    nD = max(lJ) - 1
    
    isZ = !is.null(z)
    
    if(isZ) {
      if(!is.matrix(z)) {
        z = as.matrix(z)
      }
      
      nB = ncol(z)
      
      if(type == "1p") {
        
        npar = n*nT + (nT - 1) + nD*J + nB
        
        post_names = c(paste0("lambda", 1:(nT - 1)),
                       as.vector(sapply(1:nT, function(x) {paste0(paste0("theta", x, "_"), 1:n)})),  
                       as.vector(sapply(1:J, function(x) {paste0(paste0("d", x, "_"), 1:nD)})),
                       paste0("beta", 1:nB))
        
        post = matrix(nrow = iter - burn, ncol = npar)
        post[1, ] = start_val(x = start, type = type , isZ = isZ, n = n, nT = nT, 
                              nD = nD, J = J, nB = nB)
        
      } else if(type == "2p") { 
        
        npar = n*nT + (nT - 1) + J + nD*J + nB
        
        post_names = c(paste0("lambda", 1:(nT - 1)),
                       as.vector(sapply(1:nT, function(x) {paste0(paste0("theta", x, "_"), 1:n)})),  
                       as.vector(sapply(1:J, function(x) {paste0(paste0("d", x, "_"), 1:nD)})),
                       paste0("a", 1:J),
                       paste0("beta", 1:nB))
        
        post = matrix(nrow = iter - burn, ncol = npar)
        post[1, ] = start_val(x = start, type = type , isZ = isZ, n = n, nT = nT, 
                              nD = nD, J = J, nB = nB)
      
      }
      
      
    } 
    else {
      
      if(type == "1p") {
        
        npar = n*nT + (nT - 1) + nD*J
        
        post_names = c(paste0("lambda", 1:(nT - 1)),
                       as.vector(sapply(1:nT, function(x) {paste0(paste0("theta", x, "_"), 1:n)})),  
                       as.vector(sapply(1:J, function(x) {paste0(paste0("d", x, "_"), 1:nD)})))
        
        post = matrix(nrow = iter - burn, ncol = npar)
        post[1, ] = start_val(x = start, type = type , isZ = isZ, n = n, nT = nT, 
                              nD = nD, J = J, nB = nB)
        
      } else if(type == "2p") {
        
        
        npar = n*nT + (nT - 1) + J + nD*J
        
        post_names = c(paste0("lambda", 1:(nT - 1)),
                       as.vector(sapply(1:nT, function(x) {paste0(paste0("theta", x, "_"), 1:n)})),  
                       as.vector(sapply(1:J, function(x) {paste0(paste0("d", x, "_"), 1:nD)})),
                       paste0("a", 1:J))
        
        post = matrix(nrow = iter - burn, ncol = npar)
        post[1, ] = start_val(x = start, type = type , isZ = isZ, n = n, nT = nT, 
                              nD = nD, J = J, nB = nB)
        
        
      }
    }
    
    c.ix = calc.ix(post_names, npar)
    ix = c.ix$ix
    ixe = c.ix$ixe
    
    accept = numeric(iter)
    accept[1] = 1
    
    corr_theta = matrix(data = 0, nrow = iter - burn, ncol = nT - 1)
    
    if(!is.null(z)) {
      if(type == "1p") {
        if(ntheta == 2) {
          lt1PR2D(x = x,
                  z = z,
                  iter = iter,
                  burn = burn,
                  delta = delta,
                  post = post,
                  ix = ix,
                  ixe = ixe,
                  npar = npar,
                  n = n,
                  nB = nB,
                  J = J,
                  nDmax = nD,
                  lJ = lJ,
                  nT = nT,
                  tJ = id,
                  corr_theta = corr_theta,
                  accept = accept,
                  eps = .Machine$double.eps,
                  display_progress = progress)
        } else {
          lt1PR(x = x,
                z = z,
                iter = iter,
                burn = burn,
                delta = delta,
                post = post,
                ix = ix,
                ixe = ixe,
                npar = npar,
                n = n,
                nB = nB,
                J = J,
                nDmax = nD,
                lJ = lJ,
                nT = nT,
                tJ = id,
                corr_theta = corr_theta,
                accept = accept,
                eps = .Machine$double.eps,
                display_progress = progress)
        }
      } else if(type == "2p") {
        if(ntheta == 2) {
          lt2PR2D(x = x,
                  z = z,
                  iter = iter,
                  burn = burn,
                  delta = delta,
                  post = post,
                  ix = ix,
                  ixe = ixe,
                  npar = npar,
                  n = n,
                  nB = nB,
                  J = J,
                  nDmax = nD,
                  lJ = lJ,
                  nT = nT,
                  tJ = id,
                  corr_theta = corr_theta,
                  accept = accept,
                  eps = .Machine$double.eps,
                  display_progress = progress)
        } else {
          lt2PR(x = x,
                z = z,
                iter = iter,
                burn = burn,
                delta = delta,
                post = post,
                ix = ix,
                ixe = ixe,
                npar = npar,
                n = n,
                nB = nB,
                J = J,
                nDmax = nD,
                lJ = lJ,
                nT = nT,
                tJ = id,
                corr_theta = corr_theta,
                accept = accept,
                eps = .Machine$double.eps,
                display_progress = progress)
        }
      }
    } 
    else {
      if(type == "1p") {
        if(ntheta == 2) {
          lt1PNR2D(x = x,
                   iter = iter,
                   burn = burn,
                   delta = delta,
                   post = post,
                   ix = ix,
                   ixe = ixe,
                   npar = npar,
                   n = n,
                   J = J,
                   nDmax = nD,
                   lJ = lJ,
                   nT = nT,
                   tJ = id,
                   corr_theta = corr_theta,
                   accept = accept,
                   eps = .Machine$double.eps,
                   display_progress = progress)
        } else {
          lt1PNR(x = x,
                 iter = iter,
                 burn = burn,
                 delta = delta,
                 post = post,
                 ix = ix,
                 ixe = ixe,
                 npar = npar,
                 n = n,
                 J = J,
                 nDmax = nD,
                 lJ = lJ,
                 nT = nT,
                 tJ = id,
                 corr_theta = corr_theta,
                 accept = accept,
                 eps = .Machine$double.eps,
                 display_progress = progress)
        }
      } else if(type == "2p") {
        if(ntheta == 2) {
          lt2PNR2D(x = x,
                   iter = iter,
                   burn = burn,
                   delta = delta,
                   post = post,
                   ix = ix,
                   ixe = ixe,
                   npar = npar,
                   n = n,
                   J = J,
                   nDmax = nD,
                   lJ = lJ,
                   nT = nT,
                   tJ = id,
                   corr_theta = corr_theta,
                   accept = accept,
                   eps = .Machine$double.eps,
                   display_progress = progress)
        } else {
          lt2PNR(x = x,
                 iter = iter,
                 burn = burn,
                 delta = delta,
                 post = post,
                 ix = ix,
                 ixe = ixe,
                 npar = npar,
                 n = n,
                 J = J,
                 nDmax = nD,
                 lJ = lJ,
                 nT = nT,
                 tJ = id,
                 corr_theta = corr_theta,
                 accept = accept,
                 eps = .Machine$double.eps,
                 display_progress = progress)
        }
      }
    }
    
    d_end = cumsum(rep(nD, J))
    d_start = d_end - nD + 1
    d_00 = max(lJ - 1) - (lJ - 1)
    d_null = d_end[d_00 > 0] + ix[3]
    d_null_start = d_null - d_00[d_00 > 0]
    d_null = d_null - 1
    if(length(d_null) > 0) {
      for(j in 1:length(d_null)) {
        post[, d_null_start[j]:d_null[j]] <- 0
      }
    }
    
    colnames(x) = paste0("x", 1:J)
    colnames(post) = post_names
    
    accept.rate = mean(accept)
    
    result = list(post = post, accept = accept, accept.rate = accept.rate,
                  corr_theta = corr_theta)
    
    class(result) = c("hltObj")
    return(result)
}
