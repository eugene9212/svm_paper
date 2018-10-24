# create sin and cos
class.3 <- function(n, error, t = seq(0, 1, by = 0.05), seed = 1)
{
  # 1. linear
  # 2. non-linear
  # 3. Gaussian Process
  set.seed(seed)
  K <- 3
  
  ## Create multiple y
  y <- rep(1,n)
  quotient <- n %/% K
  
  # Check the Division part
  if (n %% K != 0) message(n, " is not divisible by ",K)
  
  a <- c(1:n)
  num <- seq(0, n, quotient)
  label <- c(1:K)
  
  for (i in 1:K){
    y[(num[i]+1):num[i+1]] <- label[i]
  }

  # spatio temporal predictor
  n.t <- length(t) 
  
  # parameter
  slp <-1 ; int <- 0 # for label 1
  Sigma <- diag(x = error, nrow = n.t) # for label 3
  
  # data list
  x.list <- as.list(1:n)
  
  # 1. linear
  for (i in which(y == 1)){
    temp1 <- slp*t + int
    sample <- temp1 + rnorm(length(t), mean = 0, sd = error) # sd = 0.3
    x.list[[i]] <- sample
  }
  
  # 2. non-linear
  for (i in which(y == 2)) {
    sin.m <- sin(t)
    s.mu <- sin.m + rnorm(t, mean = 0, sd = error) # sd=0.3
    x.list[[i]] <- s.mu
  }
  
  # 3. Gaussisan process
  for (i in which(y == 3)) {
    mu.t <- z.t <- matrix(0, 1, n.t)
    mu.t <- sin(t)
    z.t <- rmvnorm(1, mu.t, Sigma)
    x.list[[i]] <- z.t
  }

  obj <- list(x = x.list, y = y, t = t)
  return(obj)
} 