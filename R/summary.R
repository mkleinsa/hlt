#' @exportS3Method summary hltObj
summary.hltObj = function(object, ...) {
  args = list(...)
  param = args$param
  dimension = args$dimension
  if("digits" %in% names(args)) {
    digits = args$digits
  } else {
    digits = 3
  }
  if("transpose" %in% names(args)) {
    transpose = args$transpose
  } else {
    transpose = TRUE
  }
  cor.theta = args$cor.theta
  post = object$post
  nms = colnames(post)
  if (param == "all") {
    smry = apply(post, 2, smy, digits = digits)
  } else if (param == "lambda") {
    smry = apply(post[, grepl("lambda", nms), drop = FALSE], 2, smy, digits = digits)
  } else if (param == "alpha") {
    smry = apply(post[, grepl("^[a]", nms), drop = FALSE], 2, smy, digits = digits)
  } else if (param == "delta") {
    smry = apply(post[, grepl("^[d]", nms), drop = FALSE], 2, smy, digits = digits)
  } else if (param == "beta") {
    smry = apply(post[, grepl("beta", nms), drop = FALSE], 2, smy, digits = digits)
  } else if (param == "cor.theta") {
    smry = apply(object$corr_theta, 2, smy, digits = digits)
  } else if (param == "theta") {
    if("dimension" %in% names(args)) {
      dimension = args$dimension
      nT = object$nT
    } else {
      warning("Since no dimension argument was specified, summaries are returned
              for the general latent dimension.")
      dimension = nT
    }
    total_theta = nrow(mod$theta)
    n_per_theta = total_theta / nT
    n_per_theta * dimension
    smry = t(object$theta[((n_per_theta * (dimension - 1)) + 1):(n_per_theta * dimension), ])
  }
  
  if(transpose == FALSE) {
    return(smry)
  } else {
    return(t(smry))
  }
}

smy = function(x, digits) {
  round(c(mean = mean(x), se = sd(x), quantile(x, probs = c(0.025, 0.5, 0.975))),
        digits = digits)
}
