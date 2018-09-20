# create linear data
linear.par <- function(n, beta, t = seq(0, 1, by = 0.05), seed = 1)
{
  # set seed
  set.seed(seed)
  
  # data storage(empty list)
  u.list <- as.list(1:n)
  
  # parameter for linear model
  slope = 3; intercept = 5
  slp <- c(slope)
  int1 <- c(intercept)
  int2 <- c(intercept + beta)
  t <- matrix(t, 1, length(t))

  # create y
  y <- rep(1, n)
  y[1:(n/2)] <- -1
  
  # create x.list
  for (i in 1:(n/2)) {

    temp1 <- t * slp + int1
    m.sample <- temp1 + rnorm(length(t), mean = 0, sd = 0.7)
    u.list[[i]] <- m.sample
  }

  for (i in (n/2+1):n) {
    
    temp2 <- t * slp + int2
    p.sample <- temp2 + rnorm(length(t), mean = 0, sd = 0.7)
    
    u.list[[i]] <- p.sample
  }
  
  x.list <- u.list
  obj <- list(x = x.list, y = y, t = t)
  return(obj)
}