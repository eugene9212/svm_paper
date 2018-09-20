# create Gaussian Process (1-dimension)
gp.1dim.wcov <- function(n, beta, t = seq(0, 1, by = 0.05), seed = 1)
{
  set.seed(seed)
  
  # binary response
  y <- rep(1, n)
  y[1:(n/2)] <- -1
  
  # spatio temporal predictor
  n.t <- length(t) 

  # covariance matrix of gaussian process
  Sigma <- 0.1 * outer(t, t, function(a, b) exp(-abs(a-b))) 
  
  u.list <- as.list(1:n)
  
  for (i in 1:n) {
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
  
  # observed covariate process, 
  # y = 1: X(t) = U(t) + beta
  # y =-1: X(t) = U(t)
  
  x.list <- u.list
  for (i in 1:n) {
    if (y[i] == 1) {
      x.list[[i]] <-  u.list[[i]] + beta
    } else x.list[[i]] <- u.list[[i]]
  }
  obj <- list(x = x.list, y = y, t = t)
  return(obj)
}
