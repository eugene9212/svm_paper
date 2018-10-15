# create Gaussian Process (1-dimension)
gp.1dim.sc.K <- function(n, K, t = seq(0, 1, by = 0.05), seed = 1)
{
  set.seed(seed)

  # create y
  y <- rep(1,n)
  quotient <- n %/% K
  
  # Check the Division part
  if (n %% K != 0) message(n, " is not divisible by ",K)
  
  num <- seq(0, n, quotient)
  label <- c(1:K)
  
  for (i in 1:K){
    y[(num[i]+1):num[i+1]] <- label[i]
  }
  
  plus <- K%/%2 
  minus <- c(K-plus)
  
  # spatio temporal predictor
  n.t <- length(t) 

  # covariance matrix of gaussian process
  Sigma <- diag(n.t)
  
  u.list <- as.list(1:n)
  
  
  # create x.list
  mu.t <- matrix(0, 1, n.t) # underlying covariate process, Z(t)
  mu.t <- sin(t) # mean vector of gaussian process
  
  for (i in 1:plus){
    for (j in (num[i]+1):num[i+1]){
      z.t <- matrix(0, 1, n.t) # underlying covariate process, Z(t)
      z.t <- rmvnorm(1, mu.t, Sigma)
      u.list[[j]] <- z.t
    }
  }
  
  mu.t <- matrix(0, 1, n.t) # underlying covariate process, Z(t)
  mu.t <- cos(t) # mean vector of gaussian process
  for (i in 1:minus){
    for (j in (num[i+plus]+1):num[i+plus+1]){
      z.t <- matrix(0, 1, n.t) # underlying covariate process, Z(t)
      z.t <- rmvnorm(1, mu.t, Sigma)
      u.list[[j]] <- z.t
      }
  }
  
  x.list <- u.list

  obj <- list(x = x.list, y = y, t = t)
  return(obj)
}

