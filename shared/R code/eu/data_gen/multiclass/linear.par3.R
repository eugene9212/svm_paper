# create linear data
linear.par3 <- function(n, beta, t = seq(0, 1, by = 0.05), seed = 1)
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
  int3 <- c(intercept + 2*beta)
  t <- matrix(t, 1, length(t))
  
  # create y
  y <- c()
  y[1:(n/3)] <- 1
  y[(n/3+1):(2*n/3)] <- 2
  y[(2*n/3+1):(n+1)] <- 3
  
  # create x.list
  for (i in 1:(n/3)) {
    temp1 <- t * slp + int1
    sample <- temp1 + rnorm(length(t), mean = 0, sd = 0.7)
    u.list[[i]] <- sample
  }
  for (i in (n/3+1):(2*n/3)) {
    temp2 <- t * slp + int2
    sample <- temp2 + rnorm(length(t), mean = 0, sd = 0.7)
    u.list[[i]] <- sample
  }
  for (i in (2*n/3+1):(n+1)) {
    temp3 <- t * slp + int3
    sample <- temp3 + rnorm(length(t), mean = 0, sd = 0.7)
    u.list[[i]] <- sample
  }
  
  x.list <- u.list
  obj <- list(x = x.list, y = y, t = t)
  return(obj)
}