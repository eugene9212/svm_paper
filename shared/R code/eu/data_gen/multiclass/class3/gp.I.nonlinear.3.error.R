# create Gaussian Process (1-dimension)
gp.I.nonlinear.3.error <- function(n, error, t = seq(0, 1, by = 0.05), seed = 1)
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
  
  # Create mean vector
  mu.t1 <- sin(3*t)
  mu.t2 <- -sin(3*t)+1
  mu.t3 <- sin(1.5*t)
  
  # x.list
  for (i in idx1) {
    # mean vector of gaussian process
    z.t <- rmvnorm(1, mu.t1, Sigma)
    x.list[[i]] <- z.t
  }
  
  for (i in idx2) {
    # mean vector of gaussian process
    z.t <- rmvnorm(1, mu.t2, Sigma)
    x.list[[i]] <- z.t
  }
  
  for (i in idx3) {
    # mean vector of gaussian process
    z.t <- rmvnorm(1, mu.t3, Sigma)
    x.list[[i]] <- z.t
  }
  
  # Calculate the True p
  true.p <- rep(0,n)
  for(i in 1:n){
    a <- dmvnorm(x=x.list[[i]], mean = mu.t1,log=TRUE)
    b <- dmvnorm(x=x.list[[i]], mean = mu.t2,log=TRUE)
    c <- dmvnorm(x=x.list[[i]], mean = mu.t3,log=TRUE)
    total <- exp(a)+exp(b)+exp(c)
    if (i %in% idx1) true.p[i] <- exp(a)/total
    else if (i %in% idx2) true.p[i] <- exp(b)/total
    else if (i %in% idx3) true.p[i] <- exp(c)/total
    else warning(paste0("No class was assigned for obs ",i))
  }
  
  obj <- list(x = x.list, y = y, t = t, true.p = true.p)
  return(obj)
}
