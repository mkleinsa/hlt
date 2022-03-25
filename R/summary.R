#' @exportS3Method summary hltObj
summary.hltObj = function(x, param, dimension = NULL, cor.theta = NULL) {
  post = x$post
  nms = colnames(post)
  if (param == "all") {
    smry = apply(post, 2, smy)
  } else if (param == "lambda") {
    smry = apply(post[, grepl("lambda", nms), drop = FALSE], 2, smy)
  } else if (param == "alpha") {
    smry = apply(post[, grepl("^[a]", nms), drop = FALSE], 2, smy)
  } else if (param == "delta") {
    smry = apply(post[, grepl("^[d]", nms), drop = FALSE], 2, smy)
  } else if (param == "beta") {
    smry = apply(post[, grepl("beta", nms), drop = FALSE], 2, smy)
  } else if (param == "cor.theta") {
    smry = apply(x$corr_theta, 2, smy)
  } else if (param == "theta") {
    smry = apply(post[, grepl(paste0("theta", dimension), nms), drop = FALSE], 2, mean)
  }
  return(smry)
}

smy = function(x) {
  c(mean = mean(x), se = sd(x), quantile(x, probs = c(0.025, 0.5, 0.7, 0.975)))
}
