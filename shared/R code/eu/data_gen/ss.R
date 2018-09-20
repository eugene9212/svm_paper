# create sin and cos
ss <- function(n, beta, t = seq(0, 1, by = 0.05), seed = 1)
{
  set.seed(seed)
  
  y <- rep(1, n)
  y[1:(n/2)] <- -1
  
  # spatio temporal predictor
  n.t <- length(t) 
  
  u.list <- as.list(1:n)
  
  for (i in 1:n) {
    
    sin.m <- sin(t)
    s.mu <- sin.m + rnorm(t, mean = 0, sd = 0.7) # sd=0.3
    
    u.list[[i]] <- s.mu
  }

  x.list <- u.list
  for (i in 1:n) {
    if (y[i] == 1) {
      x.list[[i]] <-  u.list[[i]] + beta
    } else x.list[[i]] <- u.list[[i]]
  }
  obj <- list(x = x.list, y = y, t = t)
  return(obj)
}