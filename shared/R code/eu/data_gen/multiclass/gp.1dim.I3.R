# create Gaussian Process (1-dimension)
gp.1dim.I3 <- function(n, beta, t = seq(0, 1, by = 0.05), seed = 1)
{
  set.seed(seed)
  
  # binary response
  y <- rep(1, n)
  y[(n/3+1):(2*n/3)] <- 2
  y[(2*n/3+1):(n+1)] <- 3
  
  # spatio temporal predictor
  n.t <- length(t) 
  
  # covariance matrix of gaussian process
  Sigma <- diag(n.t) 
  
  u.list <- as.list(1:n)
  
  for (i in 1:n) {
    
    mu.t <- z.t <- matrix(0, 1, n.t)
    
    # mean vector of gaussian process
    mu.t <- sin(t)
    z.t <- rmvnorm(1, mu.t, Sigma)
    
    u.list[[i]] <- z.t
  }
  
  # observed covariate process, 
  # y = 1: X(t) = U(t) + beta
  # y =-1: X(t) = U(t)
  
  x.list <- u.list
  for (i in 1:n) {
    if (y[i] == 1) {
      x.list[[i]] <-  u.list[[i]]
    } else if (y[i] == 2) {
      x.list[[i]] <-  u.list[[i]] + beta
    } else if (y[i] == 3) {
      x.list[[i]] <-  u.list[[i]] + 2*beta
    }
  }
  
  
  obj <- list(x = x.list, y = y, t = t)
  return(obj)
}
