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
  if(missing(y)) {
    post = do.call(rbind, Map(f = function(x) {x$post}, mod))
    return(post)
  } else {
    
  }
}
