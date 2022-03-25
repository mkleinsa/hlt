# 
# \dontrun{
# 
# set.seed(7794)
# xdat = hltsim(n = 150, ntheta = ntheta, lambda = lambda, 
#               tJ = id, dL = dL, mua = 1, mud = 1.4, siga = 0.8, sigd = 1,
#               beta = 1)
# x = xdat$x
# mod = hlt(x, ntheta = 3, id = id, iter = 40000, burn = 30000, delta = 0.05)
# mod$accept.rate
# post = mod$post
# 
# # Regression example 
# set.seed(2034)
# nB = 1
# n = 300
# z = matrix(sample(0:1, n, replace = TRUE), nrow = n, ncol = nB)
# xdat = hltsim(n = n, ntheta = 3, lambda = c(0.8, 0.4, 0.3), 
#               tJ = id, dL = 4, mua = 1, mud = 1.4, siga = 0.8, sigd = 1,
#               regression = TRUE, z = z, nB = nB, beta = 1)
# apply(xdat$x, 2, table)
# lm(xdat$theta[,4] ~ xdat$z)
# xdat$s.beta
# x = xdat$x
# x[,1] = ifelse(x[,1] == 3, 2, x[,1])
# x[,1] = ifelse(x[,1] == 2, 1, x[,1])
# x[,4] = ifelse(x[,4] == 3, 2, x[,4])
# x[,5] = ifelse(x[,5] == 3, 2, x[,5])
# mod = hlt(x, z = z, id = id, iter = 1e5, delta = 0.017)
# mod = hlt(x, z = z, id = id, iter = 1e6, burn = 9e5, delta = 0.017)
# 
# mod = hlt(x, z = z, id = id, iter = 2.2e6, burn = 2e6, delta = 0.013)
# mod$accept.rate
# post = mod$post
# apply(post, 2, mean)
# 
# smy = function(x) {c(mean = mean(x), se = sd(x), quantile(x, probs = c(0.025, 0.5, 0.7, 0.975)))}
# apply(post[, "beta1", drop = FALSE], 2, smy)
# 
# summary(mod, param = "all")
# summary(mod, param = "beta")
# summary(mod, param = "lambda")
# summary(mod, param = "alpha")
# summary(mod, param = "delta")
# summary(mod, param = "theta", dimension = 1)
# summary(mod, param = "theta", dimension = 2)
# summary(mod, param = "theta", dimension = 3)
# summary(mod, param = "theta", dimension = 4)
# 
# cor(as.vector(xdat$s.delta[-1,]), 
#     as.vector(matrix(summary(mod, param = "delta")[1,], nrow = 3)))
#     
# cor(summary(mod, param = "alpha")[1,], xdat$s.alpha)
# 
# summary(mod, param = "cor.theta", cor.theta = c(1,2))
# xdat$s.lam.cor
# 
# th = summary(mod, param = "theta", dimension = 1)
# cor(th, xdat$theta[,1])
# 
# plot(mod, "lambda1", type = "trace")
# plot(mod, "lambda2", type = "trace")
# plot(mod, "lambda3", type = "trace")
# plot(mod, "a1", type = "trace")
# plot(mod, "d2", type = "trace")
# plot(mod, "beta1", type = "trace")
# 
# plot(mod, 1, type = "icc")
# plot(mod, 2, type = "icc")
# plot(mod, 3, type = "icc")
# plot(mod, 4, type = "icc")
# plot(mod, 5, type = "icc")
# plot(mod, 6, type = "icc")
# plot(mod, 7, type = "icc", min = -10, max = 10)
# 
# # example 
# data("asti")
# x = as.matrix(asti[, 1:25]) - 1
# z = asti[, 26:27]
# z[, 1] = (z[, 1] == "students") * 1
# z[, 2] = (z[, 2] == "male") * 1
# z = as.matrix(z)
# id = c(0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4)
# mod = hlt(x, z = z, id = id, iter = 2.2e6, burn = 2e6, delta = 0.013)
# 
# 
# ### Test each combination of input values ###
# n = 200
# nB = 5
# z = matrix(sample(0:1, n, replace = TRUE), nrow = n, ncol = nB)
# 
# # PCM, No regression, 2 dimensions
# xdat = hltsim(n = n, ntheta = 2, lambda = c(0.7, 0.7), 
#               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
#               dL = 2, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#               regression = FALSE, z = NULL)
# apply(xdat$x, 2, table)
# mod = hlt(xdat$x, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1), 
#           iter = 200, burn = 100, delta = 0.013)
# 
# # PCM, Regression, 2 dimensions
# xdat = hltsim(n = n, ntheta = 2, lambda = c(0.7, 0.7), 
#               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
#               dL = 2, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#               regression = FALSE, z = z)
# apply(xdat$x, 2, table)
# mod = hlt(xdat$x, z = z, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1), 
#           iter = 200, burn = 100, delta = 0.013)
# 
# # PCM, No regression, 4 dimensions
# xdat = hltsim(n = n, ntheta = 4, lambda = c(0.7, 0.2, 1.2, 0.4), 
#               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#                      3,3,3,3,3,3),
#               dL = 2, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#               regression = FALSE, z = NULL)
# apply(xdat$x, 2, table)
# mod = hlt(xdat$x, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#           3,3,3,3,3,3), iter = 200, burn = 100, delta = 0.013)
# 
# # PCM, Regression, 4 dimensions
# xdat = hltsim(n = n, ntheta = 4, lambda = c(0.7, 0.2, 1.2, 0.4), 
#               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#                      3,3,3,3,3,3),
#               dL = 2, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#               regression = FALSE, z = z)
# apply(xdat$x, 2, table)
# mod = hlt(xdat$x, z = z, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#           3,3,3,3,3,3), iter = 200, burn = 100, delta = 0.013)
# 
# 
# 
# 
# 
# 
# # GPCM, No regression, 2 dimensions
# xdat = hltsim(n = n, ntheta = 2, lambda = c(0.7, 0.7), 
#               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
#               dL = 10, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#               regression = FALSE, z = NULL)
# apply(xdat$x, 2, table)
# mod = hlt(xdat$x, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1), 
#           iter = 200, burn = 100, delta = 0.013)
# 
# # GPCM, Regression, 2 dimensions
# xdat = hltsim(n = n, ntheta = 2, lambda = c(0.7, 0.7), 
#               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),
#               dL = 10, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#               regression = FALSE, z = z)
# apply(xdat$x, 2, table)
# mod = hlt(xdat$x, z = z, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1), 
#           iter = 200, burn = 100, delta = 0.013)
# 
# # GPCM, No regression, 4 dimensions
# xdat = hltsim(n = n, ntheta = 4, lambda = c(0.7, 0.2, 1.2, 0.4), 
#               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#                      3,3,3,3,3,3),
#               dL = 10, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#               regression = FALSE, z = NULL)
# apply(xdat$x, 2, table)
# mod = hlt(xdat$x, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#           3,3,3,3,3,3), iter = 200, burn = 100, delta = 0.013)
# 
# # GPCM, Regression, 4 dimensions
# xdat = hltsim(n = n, ntheta = 4, lambda = c(0.7, 0.2, 1.2, 0.4), 
#               tJ = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#                      3,3,3,3,3,3),
#               dL = 10, mua = NULL, mud = 1.4, siga = NULL, sigd = 1,
#               regression = FALSE, z = z)
# apply(xdat$x, 2, table)
# mod = hlt(xdat$x, z = z, id = c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2, 
#           3,3,3,3,3,3), iter = 200, burn = 100, delta = 0.013)
# 
# }
# 