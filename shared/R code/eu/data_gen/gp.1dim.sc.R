# create Gaussian Process (1-dimension)
gp.1dim.sc <- function(n, t = seq(0, 1, by = 0.05), seed = 1)
{
  set.seed(seed)
  # binary response
  y <- rep(1, n)
  y[1:(n/2)] <- -1
  
  # spatio temporal predictor
  n.t <- length(t) 

  # covariance matrix of gaussian process
  Sigma <- diag(n.t)
  
  u.list <- as.list(1:n)
  for (i in 1:(n/2)) {
    mu.t <- z.t <- matrix(0, 1, n.t)
    # underlying covariate process, Z(t)
    
    # mean vector of gaussian process
    mu.t <- sin(t)
    z.t <- rmvnorm(1, mu.t, Sigma)
    
    # spatially correlated process, U(t) = sum_neighbor Z(t)
    u.t <- matrix(0, 1, n.t)
    u.t <- apply(z.t, 2, mean)
    u.list[[i]] <- u.t
  }
  
  for (i in (n/2+1):n) {
    mu.t <- z.t <- matrix(0, 1, n.t)
    # underlying covariate process, Z(t)
    
    # mean vector of gaussian process
    mu.t <- cos(t)
    z.t <- rmvnorm(1, mu.t, Sigma)
    
    # spatially correlated process, U(t) = sum_neighbor Z(t)
    u.t <- matrix(0, 1, n.t)
    u.t <- apply(z.t, 2, mean)
    u.list[[i]] <- u.t
  }

  x.list <- u.list

  obj <- list(x = x.list, y = y, t = t)
  return(obj)
}

