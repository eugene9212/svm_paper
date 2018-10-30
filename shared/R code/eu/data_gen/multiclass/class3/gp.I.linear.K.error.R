# create Gaussian Process (1-dimension)
gp.I.linear.K.error <- function(n, error, beta, K, t = seq(0, 1, by = 0.05), seed = 1)
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
  
  # covariance matrix of gaussian process
  Sigma <- error*diag(n.t) 
  
  x.list <- u.list <- as.list(1:n)
  
  for (i in 1:n) {
    
    mu.t <- z.t <- matrix(0, 1, n.t)
    
    # mean vector of gaussian process
    mu.t <- t
    z.t <- rmvnorm(1, mu.t, Sigma)
    
    u.list[[i]] <- z.t
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
