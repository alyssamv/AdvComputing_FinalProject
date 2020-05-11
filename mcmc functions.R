
ll_f = function(X, Y, beta, Sigma, rho) {
  sig.vec = rep(Sigma, 7)
  invsigma = diag(1/sig.vec)
  mu_X = (X) %*% diag(beta) + rho*Y
  
  a = log(Sigma)
  b = (X - (mu_X)) %*% invsigma %*% t(X - (mu_X)) %>% as.numeric
  c = 7*log(2*pi)
  
  l = -0.5*(as.vector(a + b) + c)
  return(l)
}


ll_rho = function(rhos) {
  # Prior for correlation rho between time points
  prob = dtruncnorm(rhos, a = 0, b = 1, mean = 0.5, sd = sqrt(0.2))
  return(log(prob))
}


ll_invsigma = function(sigma) {
  # Prior for the inverse of the error covariance matrix
  return(log(dinvgamma(1/sigma, shape = 0.0001, rate = 0.0001)))
}


ll_beta = function(betas, sigma) {
  return((log(dnorm(betas, mean = 0, sd = sqrt(sigma)))))
}



logpost_all = function(X, Y, beta, rho, Sigma){  
  post_beta = ll_beta(beta, Sigma) + ll_invsigma(Sigma) + ll_f(X, Y, beta, Sigma, rho)
  post_rho = ll_rho(rho) + ll_f(X, Y, beta, Sigma, rho)
  post_sigma = ll_invsigma(Sigma) + ll_f(X, Y, beta, Sigma, rho)
  
  return(c(post_beta, post_rho, post_sigma))
}



MH_step = function(train_X, train_Y, params, b_avec, r_a, sig_a){
  
  # cast X as matrix 
  train_X = unname(train_X) %>% as.matrix
  
  # initialize vectors for posterior estimates
  bb = c() # beta vec
  rr = NA # rho
  ss = NA # Sigma
  
  # current values 
  p_betas = cur_betas = params[1:7]
  cur_rho = params[8]
  cur_sigma = params[9]
  
  ### Component-wise MH 
  # estimate Betas
  bnum = vector(length = 7)
  bdenom = vector(length = 7)
  for (j in 1:7) {
    cur_betas[j] = p_betas[j]
    p_betas[j] = cur_betas[j] + (runif(1) - 0.5) * 2 * b_avec[j] # use random walk to sample next value
    
    bnum = logpost_all(train_X, train_Y, beta = p_betas, rho = cur_rho, Sigma = cur_sigma)
    bdenom = logpost_all(train_X, train_Y, beta = cur_betas, rho = cur_rho, Sigma = cur_sigma)
    if (log(runif(1)) < (bnum[j] - bdenom[j])) {
      bb[j] = p_betas[j]
    } else {
      bb[j] = cur_betas[j]}
  }
  
  #sampling rhos
  p_rho = cur_rho + (runif(1) - 0.5) * 2 * r_a
  rnum = logpost_all(train_X, train_Y, beta = bb, rho = p_rho, Sigma = cur_sigma)
  rdenom = logpost_all(train_X, train_Y, beta = bb, rho = cur_rho, Sigma = cur_sigma)
  #cur_rho = p_rho
  if (log(runif(1)) < (rnum[8] - rdenom[8])) {
    rr = p_rho
  } else {
    rr = cur_rho
  }
  
  #sampling sigmas
  p_sigma = cur_sigma + (runif(1) - 0.5) * 2 * sig_a
  
  # if (p_sigma < 0) {
  #   snum = Inf
  # } else {
  #   snum = logpost_all(train_X, train_Y, beta = bb, rho = rr, Sigma = p_sigma)
  # }
  snum = logpost_all(train_X, train_Y, beta = bb, rho = rr, Sigma = p_sigma)
  sdenom = logpost_all(train_X, train_Y, beta = bb, rho = rr, Sigma = cur_sigma)
  if (log(runif(1)) < (snum[9] - sdenom[9])) {
    ss = p_sigma
  } else {
    ss = cur_sigma}
  
  return(c(bb, rr, ss))
}


numunique <- function(mat){
  apply(mat, 2, function(i){
    length(unique(i))
  })
}