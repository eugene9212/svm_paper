# create Gaussian Process (1-dimension)
gp.1dim.AR <- function(n, beta, p = 5, rho = 0.5, t = seq(0, 1, by = 0.05), seed = 1)
{
  set.seed(seed)
  
  # binary response
  y <- rep(1, n)
  y[1:(n/2)] <- -1
  
  # spatio temporal predictor
  n.t <- length(t) 
  
  # create AR(p) covariance structure
  rho.1 <- rho^c(0:p)
  n.t <- length(t)
  A <- matrix(c(rho.1, rep(0, n.t-p)), n.t, n.t, byrow = T)
  makeSymm <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    return(m)
  }
  Sigma <- makeSymm(A)
  
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
      x.list[[i]] <-  u.list[[i]] + beta
    } else x.list[[i]] <- u.list[[i]]
  }
  obj <- list(x = x.list, y = y, t = t)
  return(obj)
}
