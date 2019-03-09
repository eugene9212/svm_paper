linear.cross.K <- function(n, K, t = seq(0, 1, by = 0.05), seed = 1)
{
  # set seed
  set.seed(seed)
  
  # data storage(empty list)
  u.list <- as.list(1:n)
  
  # parameter for linear model
  plus <- K%/%2 
  minus <- c(K-plus)
  slp1 <- runif(plus, min = 0, max = 1.5)
  slp2 <- -runif(minus, min = 0, max = 1.5)
  int1 <- runif(plus, min = 4, max = 6)
  int2 <- runif(minus, min = 4, max = 6)
  t <- matrix(t, 1, length(t))
  
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

  # create x.list
  for (i in 1:plus){
    for (j in (num[i]+1):num[i+1]){
      temp1 <- slp1[i]*t + int1[i]
      sample <- temp1 + rnorm(length(t), mean = 0, sd = 0.7) # sd = 0.3
      u.list[[j]] <- sample
    }
  }
    
  for (i in 1:minus){
    for (j in (num[i+plus]+1):num[i+plus+1]){
      print(j)
      temp2 <- slp2[i]*t + int2[i]
      sample <- temp2 + rnorm(length(t), mean = 0, sd = 0.7) # sd = 0.3
      u.list[[j]] <- sample
    }
  }
    
  x.list <- u.list
  
  obj <- list(x = x.list, y = y, t = t)
  return(obj)
}