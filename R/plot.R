#' 
#' @import ggplot2
#' @exportS3Method plot hltObj
plot.hltObj = function(mod, type, x, ...) {
  post = mod$post
  nr = nrow(post)
  if(type == "trace") {
    ggplot(data.frame(x = 1:nr, y = post[, x]), aes(x, y)) + geom_line() + 
      xlab("iteration") + ylab("value") + get_theme()
  } else if(type == "icc") {
    plot.hltObj.icc(mod, x = x, type = type, ...)
  }
}

get_theme = function() {
  theme_bw() + 
  theme(panel.grid.major = element_line(colour="#DDDDDD", size = (.5)),
        panel.grid.minor = element_line(size = (0.2), colour="#DDDDDD"),
        panel.border = element_rect(colour = "black", size = 1.5),
        axis.ticks = element_line(size = 1), axis.ticks.length = unit(.2, "cm"),
        plot.title = element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.key.size = unit(1, 'lines'),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 15, margin = margin(t = 12, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1))
}

plot.hltObj.icc = function(mod, x, type, ...) {
  args = list(...)
  if(type == "icc") {
    plt = icc_curve(mod, x, min = ifelse(is.null(args$min), -4, args$min),
                    max = ifelse(is.null(args$max), 4, args$max))
  }
  return(plt)
}

#' @importFrom tidyr pivot_longer
icc_curve = function(mod, x, min = -4, max = 4) {
  if(any(grepl("^[a]", colnames(mod$post)))) {
    alpha = as.vector(summary(mod, param = "alpha")["mean",])
  }
  kappa = summary(mod, param = "delta")["mean",]
  kappa_names = names(kappa)
  kappa_id = as.numeric(gsub(".*[d]([^.]+)[_].*", "\\1", kappa_names))
  kappa_list = vector(length = length(unique(kappa_id)), mode = "list")
  for(i in 1:length(kappa_list)) {
    kappa_list[[i]] = kappa[kappa_id == i]
  }
  pxk = function(theta) {
    stat = "mean"
    return(as.numeric(sapply(theta, function(theta) {
      1 / (1 + exp(-alpha[x] * (theta - kappa_list[[x]])))
    })))
  }
  if(length(kappa_list[[x]]) == 1) {
    pdata = data.frame(x = seq(min, max, by = 0.1), y = sapply(seq(min, max, by = 0.1), pxk))
  } else {
    pdata = data.frame(x = seq(min, max, by = 0.1), y = t(sapply(seq(min, max, by = 0.1), pxk)))
  }
  kappa_nms = names(kappa_list[[x]])
  names(pdata) = c("x", kappa_nms)
  pdata = tidyr::pivot_longer(pdata, cols = kappa_nms)
  plt = ggplot(pdata, aes(x, value, color = name)) + get_theme() + 
    xlab("theta") + ylab("p") + geom_line()
  return(plt)
}
