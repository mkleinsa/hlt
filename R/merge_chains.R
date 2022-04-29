#' Merge Chains from hlt method
#'
#' @param x object of class "hltObjList"
#' @param ... other arguments
#' @export
merge_chains = function (x, ...) {
  UseMethod("merge_chains", x)
}

#' @exportS3Method merge_chains hltObjList
merge_chains.hltObjList = function(x, ...) {
  if(missing(...)) {
    nchains = length(x)
    post = do.call(rbind, Map(f = function(y) {y$post}, x))
    nT = x[[1]]$nT
    thetas = Map(f = function(y) {y$theta}, x)
    means = as.data.frame(Map(f = function(y) {y$mean}, thetas))
    sds = as.data.frame(Map(f = function(y) {y$sd}, thetas))
    rmeans = rowMeans(means)
    rsts = rowMeans(sds) 
    theta = data.frame(mean = rmeans, sd = rsts)
    result = list(post = post, 
                  theta = theta,
                  nT = nT,
                  nchains = nchains)
    class(result) = c("hltObj")
    return(result)
  } else {
    
  }
}
