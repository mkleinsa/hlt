
auto_start = function (x) {
  list(
    lambda = summary(x, param = "lambda")[, "mean"] + rnorm(5, 0, 0.1),
    theta = x$theta[, "mean"] + rnorm(6774, 0, 0.1), 
    delta = summary(x, param = "delta")[, "mean"] + rnorm(75, 0, 0.1), 
    alpha = summary(x, param = "alpha")[, "mean"] + abs(rnorm(25, 0, 0.1)), 
    beta = c()
  )
}