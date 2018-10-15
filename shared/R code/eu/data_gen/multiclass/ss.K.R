# create sin and cos
ss.K <- function(n, beta, K, t = seq(0, 1, by = 0.05), seed = 1)
{
  set.seed(seed)
  
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
  
  beta1 <- seq(0, (K-1), 1) * beta
  
  # spatio temporal predictor
  n.t <- length(t) 
  
  u.list <- as.list(1:n)
  
  for (i in 1:n) {
    sin.m <- sin(t)
    s.mu <- sin.m + rnorm(t, mean = 0, sd = 0.7) # sd=0.3
    
    u.list[[i]] <- s.mu
  }
  
  ## Create x.list
  for (i in 1:K){
    for (j in (num[i]+1):num[i+1]){
      x.list[[j]] <-  u.list[[j]] + beta1[i]
    }
  }
  
  obj <- list(x = x.list, y = y, t = t)
  return(obj)
}