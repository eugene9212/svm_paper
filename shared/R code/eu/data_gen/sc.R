# create sin and cos
sc <- function(n, t = seq(0, 1, by = 0.05), seed = 1)
{
  set.seed(seed)
  
  u.list <- as.list(1:n)
  n.t <- length(t) 

  # create y
  y <- rep(1, n)
  y[1:(n/2)] <- -1
 
  # create x
  for (i in 1:(n/2)) {
    sin.m <- sin(t)
    s.mu <- sin.m + rnorm(t, mean = 0, sd = 0.7)

    u.list[[i]] <- s.mu
  }
  
  for (i in (n/2+1):n) {
    cos.m <- cos(t)
    c.mu <- cos.m + rnorm(t, mean = -0.5, sd = 0.7)
    
    u.list[[i]] <- c.mu
  }
  
  x.list <- u.list
  
  obj <- list(x = x.list, y = y, t = t)
  return(obj)
}