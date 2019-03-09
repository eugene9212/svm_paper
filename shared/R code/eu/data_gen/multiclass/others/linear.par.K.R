# create linear data
linear.par.K <- function(n, beta, K, t = seq(0, 1, by = 0.05), seed = 1)
{
  # set seed
  set.seed(seed)
  
  # data storage(empty list)
  u.list <- as.list(1:n)
  
  # parameter for linear model
  slp <- 3
  int <- 5 + seq(0, (K-1), 1) * beta
  
  t <- matrix(t, 1, length(t))
  
  # create y
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
  
  # create x.list
  for (i in 1:K){
    # print(i)
    for (j in (num[i]+1):num[i+1]){
      temp <- t * slp + int[i]
      sample <- temp + rnorm(length(t), mean = 0, sd = 0.7)
      u.list[[j]] <- sample
    }
  }
  
  x.list <- u.list
  obj <- list(x = x.list, y = y, t = t)
  return(obj)
}
