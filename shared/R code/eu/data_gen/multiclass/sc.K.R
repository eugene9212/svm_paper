# create sin and cos
sc.K <- function(n, K, t = seq(0, 1, by = 0.05), seed = 1)
{
  set.seed(seed)
  
  u.list <- as.list(1:n)
  n.t <- length(t) 
  
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
  
  # assign sin and cos
  plus <- K%/%2 
  minus <- c(K-plus)
  
  # create x.list
  for (i in 1:plus){
    for (j in (num[i]+1):num[i+1]){
      sin.m <- sin(t)
      s.mu <- sin.m + rnorm(t, mean = 0, sd = 0.7)
      
      u.list[[j]] <- s.mu
    }
  }
  
  for (i in 1:minus){
    for (j in (num[i+plus]+1):num[i+plus+1]){
      cos.m <- cos(t)
      c.mu <- cos.m + rnorm(t, mean = -0.5, sd = 0.7)
      
      u.list[[j]] <- c.mu
    }
  }
  
  x.list <- u.list
  
  obj <- list(x = x.list, y = y, t = t)
  return(obj)
}