# create Gaussian Process (1-dimension)
gp.1dim.sc <- function(n, error, cov, t = seq(0, 1, by = 0.05), seed = 1)
{
  set.seed(seed)
  # binary response
  y <- rep(1, n)
  y[1:(n/2)] <- -1
  
  # spatio temporal predictor
  n.t <- length(t) 

  # covariance matrix of gaussian process
  if (cov=="I") {
    Sigma <- error*diag(n.t)
  } else if(cov=="AR") {
    rho <- 0.7
    idx.t <- c(1:n.t)
    Sigma <- error * outer(idx.t, idx.t, function(a, b) rho^abs(a-b))
  } else if(cov=="CS") {
    rho <- 0.3
    tmp <- matrix(rho,n.t,n.t)
    diag(tmp) <- 1
    Sigma <- error * tmp
  } else warning("covariance structure was not specified")
  
  # Divide index
  idx1 <- which(y==-1)
  idx2 <- which(y==1)
  
  # Create mean vector
  mu.t1 <- sin(t)
  mu.t2 <- cos(t)
  
  # x.list
  x.list <- as.list(1:n)
  for (i in idx1) {
    z.t <- rmvnorm(1, mu.t1, Sigma)
    x.list[[i]] <- z.t
  }
  
  for (i in idx2) {
    z.t <- rmvnorm(1, mu.t2, Sigma)
    x.list[[i]] <- z.t
  }
  
  # Calculate the True p
  true.p <- rep(0,n)
  for(i in 1:n){
    a <- dmvnorm(x=x.list[[i]], mean = mu.t1,log=TRUE)
    b <- dmvnorm(x=x.list[[i]], mean = mu.t2,log=TRUE)
    total <- exp(a)+exp(b)
    if (i %in% idx1) true.p[i] <- exp(a)/total
    else if (i %in% idx2) true.p[i] <- exp(b)/total
    else warning(paste0("No class was assigned for obs ",i))
  }
  
  obj <- list(x = x.list, y = y, t = t, true.p = true.p)
  return(obj)
}

