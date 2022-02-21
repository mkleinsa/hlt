#' @exportS3Method summary hltObj
summary.hltObj = function(x, param, dimension = NULL, cor.theta = NULL) {
  nms = colnames(post)
  if (param == "lambda") {
    smry = apply(post[, grepl("lambda", nms), drop = FALSE], 2, smy)
  } else if (param == "alpha") {
    smry = apply(post[, grepl("^[a]\\d+$", nms), drop = FALSE], 2, smy)
  } else if (param == "delta") {
    smry = apply(post[, grepl("^[d]\\d+$", nms), drop = FALSE], 2, smy)
  } else if (param == "beta") {
    smry = apply(post[, grepl("beta", nms), drop = FALSE], 2, smy)
  } else if (param == "cor.theta") {
    post[, grepl(paste0("theta", cor.theta[1]), nms)]
    post[, grepl(paste0("theta", cor.theta[2]), nms)]
    cor.post = matrix(nrow = nrow(post), ncol = 1)
    for(i in 1:nrow(cor.post)) {
      cor.post[i,1] = cor(x = post[i, grepl(paste0("theta", cor.theta[1]), nms)],
                          y = post[i, grepl(paste0("theta", cor.theta[2]), nms)])
    }
    smry = apply(cor.post, 2, smy)
  } else if (param == "theta") {
    smry = apply(post[, grepl(paste0("theta", dimension), nms), drop = FALSE], 2, mean)
  }
  
  return(smry)
}

cor.hltObj = function(x, y) {
  
}

smy = function(x) {
  c(mean = mean(x), se = sd(x), quantile(x, probs = c(0.025, 0.5, 0.7, 0.975)))
}
