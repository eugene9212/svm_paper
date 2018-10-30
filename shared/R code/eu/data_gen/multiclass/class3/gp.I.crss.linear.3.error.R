# create Gaussian Process (1-dimension)
gp.I.crss.linear.3.error <- function(n, error, t = seq(0, 1, by = 0.05), seed = 1)
{
  set.seed(seed)
  K<-3
  
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
  
  # covariance matrix of gaussian process
  Sigma <- error*diag(n.t) 
  
  x.list <- as.list(1:n)
  
  idx1 <- which(y==1)
  idx2 <- which(y==2)
  idx3 <- which(y==3)
  
  for (i in idx1) {
    mu.t <- z.t <- matrix(0, 1, n.t)
    
    # mean vector of gaussian process
    mu.t <- 3*t+1
    z.t <- rmvnorm(1, mu.t, Sigma)
    
    x.list[[i]] <- z.t
  }
  
  for (i in idx2) {
    mu.t <- z.t <- matrix(0, 1, n.t)
    
    # mean vector of gaussian process
    mu.t <- -3*t+4
    z.t <- rmvnorm(1, mu.t, Sigma)
    
    x.list[[i]] <- z.t
  }
  
  for (i in idx3) {
    mu.t <- z.t <- matrix(0, 1, n.t)
    
    # mean vector of gaussian process
    mu.t <- rep(3,n.t)
    z.t <- rmvnorm(1, mu.t, Sigma)
    
    x.list[[i]] <- z.t
  }
  
  obj <- list(x = x.list, y = y, t = t)
  return(obj)
}
