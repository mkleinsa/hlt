
get_scenario = function(isZ, type, n, nT, nD, J, nB, start) {
  
  npar_theta = n*nT
  
  theta_mean = matrix(iter - burn, nT)
  theta_mean_sq = matrix(iter - burn, nT)
  
  if(isZ) {
    if(type == "1p") {
      npar = n*nT + (nT - 1) + nD*J + nB
      npar_total = npar + npar_theta
      
      post_names = c(paste0("lambda", 1:(nT - 1)),
                     as.vector(sapply(1:nT, function(x) {paste0(paste0("theta", x, "_"), 1:n)})),  
                     as.vector(sapply(1:J, function(x) {paste0(paste0("d", x, "_"), 1:nD)})),
                     paste0("beta", 1:nB))
      
      post = matrix(nrow = iter - burn, ncol = npar)
      post[1, ] = start_val(x = start, type = type , isZ = isZ, n = n, nT = nT, 
                            nD = nD, J = J, nB = nB)
      
    } else if(type == "2p") {
      #npar_theta = n*nT
      npar = (nT - 1) + J + nD*J + nB
      npar_total = npar + npar_theta
      
      draw = matrix(nrow = 1, ncol = npar)
      draw_theta = matrix(nrow = 1, ncol = npar_theta)
      
      #theta_mean = matrix(iter - burn, nT)
      #theta_mean_sq = matrix(iter - burn, nT)
      
      post_names = c(paste0("lambda", 1:(nT - 1)),
                     as.vector(sapply(1:nT, function(x) {paste0(paste0("theta", x, "_"), 1:n)})),  
                     as.vector(sapply(1:J, function(x) {paste0(paste0("d", x, "_"), 1:nD)})),
                     paste0("a", 1:J),
                     paste0("beta", 1:nB))
      
      post_names_theta = post_names[grepl("theta", post_names)]
      post_names_no_theta = post_names[!grepl("theta", post_names)]
      
      post = matrix(nrow = iter - burn, ncol = npar)
      colnames(post) = post_names_no_theta
      
      start_values = start_val(x = start, type = type, isZ = isZ, n = n, nT = nT, 
                               nD = nD, J = J, nB = nB)
      
      draw[1, ] = start_values[post_names %in% post_names_no_theta]
      draw_theta[1, ] = start_values[post_names %in% post_names_theta]
      
      return(list(npar = npar, npar_theta = npar_theta, npar_total = npar_total, 
             draw = draw, draw_theta = draw_theta, theta_mean = theta_mean,
             theta_mean_sq = theta_mean_sq, post_names = post_names, 
             post_names_theta = post_names_theta, post_names_no_theta = post_names_no_theta,
             post = post, start_values = start_values))
    }
  } else {
    if(type == "1p") {
      npar = n*nT + (nT - 1) + nD*J
      npar_total = npar + npar_theta
      
      post_names = c(paste0("lambda", 1:(nT - 1)),
                     as.vector(sapply(1:nT, function(x) {paste0(paste0("theta", x, "_"), 1:n)})),  
                     as.vector(sapply(1:J, function(x) {paste0(paste0("d", x, "_"), 1:nD)})))
      
      post = matrix(nrow = iter - burn, ncol = npar)
      post[1, ] = start_val(x = start, type = type , isZ = isZ, n = n, nT = nT, 
                            nD = nD, J = J, nB = nB)
    } else if(type == "2p") {
      # npar = n*nT + (nT - 1) + J + nD*J
      # npar_keep = npar - n*nT
      # draw = matrix(nrow = 1, ncol = npar)
      # theta_mean = matrix(iter - burn, nT)
      # theta_mean_sq = matrix(iter - burn, nT)
      # 
      # post_names = c(paste0("lambda", 1:(nT - 1)),
      #                #as.vector(sapply(1:nT, function(x) {paste0(paste0("theta", x, "_"), 1:n)})),  
      #                as.vector(sapply(1:J, function(x) {paste0(paste0("d", x, "_"), 1:nD)})),
      #                paste0("a", 1:J))
      # 
      # start_values = start_val(x = start, type = type, isZ = isZ, n = n, nT = nT, 
      #                          nD = nD, J = J, nB = nB)
      # post = matrix(nrow = iter - burn, ncol = npar)
      # post[1, ] = start_val(x = start, type = type, isZ = isZ, n = n, nT = nT, 
      #                       nD = nD, J = J, nB = nB)
    }
  }
}