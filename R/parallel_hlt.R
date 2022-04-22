# 
# library(parallel)
# library(MASS)
# 
# starts <- rep(100, 40)
# fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
# numCores <- detectCores()
# numCores
# 
# system.time(
#   results <- lapply(starts, fx)
# )
# 
# system.time(
#   results <- mclapply(starts, fx, mc.cores = numCores)
# )
# 
# library(foreach)
# library(doParallel)
# 
# registerDoParallel(numCores)  # use multicore, set to the number of our cores
# 
# # set seed
# set.seed(153)
# 
# # load the asti data set
# data("asti")
# 
# # shift responses to range from 0 instead of 1
# x = as.matrix(asti[, 1:25]) - 1
# 
# # specify which items from x belong to each domain
# id = c(0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4)
# 
# 
# nchains = 4
# models = foreach (i = 1:nchains, .verbose = TRUE) %dopar% {
#   hlt(x, id = id, iter = 20, burn = 10, delta = 0.01, progress = TRUE)
# }
# names(models) = paste0("chain", 1:nchains)











