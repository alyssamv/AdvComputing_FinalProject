##### Log Likelihood functions

## Sigma is a vector of beta_k variances
## X is vector of observations
## mu_X??
## rho is scalar

ll_f = function(Sigma, X, mu_X) {
  invsigma = diag(1/Sigma)
  
  a = log(Sigma)
  b = (X - mu_X) %*% invsigma %*% (X - mu_X)
  c = 7*log(2*pi)
  
  l = -0.5(a + b + c)
  return(l)
}


ll_beta = function(Sigma, X) {
  invsigma = diag(1/Sigma)
  
  a = log(Sigma)
  b = X %*% invsigma %*% X
  c = 7*log(2*pi)
  
  l = -0.5(a + b + c)
  return(l)
}

ll_rho = function(rho) {
  a = log(sqrt(5))
  b = log(dnorm(sqrt(5)*(rho - 0.5)))
  c = log(pnorm(sqrt(5)/2) - pnorm(-sqrt(5)/2))
  l = a + b - c
  return(l)
}

ll_invSigma = function(Sigma) {
  invsigma = diag(1/Sigma)

  a = 0.001*log(0.001)
  b = log(factorial(0.001 - 1))
  c = 1.001*log(invsigma)
  d = 0.001/invsigma
  
  l = a - b + c + d
  return(l)
}