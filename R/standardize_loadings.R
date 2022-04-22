
standardize_loadings = function(x) {
  (x - min(x)) / (max(x) - min(x))
}