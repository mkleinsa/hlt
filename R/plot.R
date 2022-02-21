#' 
#' @import ggplot2
#' @exportS3Method plot hltObj
plot.hltObj = function(mod, x) {
  post = mod$post
  nr = nrow(post)
  ggplot(data.frame(x = 1:nr, y = post[, x]), aes(x, y)) + geom_line()
}