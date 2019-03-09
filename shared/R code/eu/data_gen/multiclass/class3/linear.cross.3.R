# create Gaussian Process
linear.cross.3 <- function(n, error, cov, rho, t = seq(0, 1, by = 0.05), seed = 1)
{
  set.seed(seed)
  K <- 3
  
  # Check the Division part
  if (n %% K != 0) message(n, " is not divisible by ",K)
  
  ##== Class labeling at y ==##
  y <- rep(1,n)
  quotient <- n %/% K
  
  a <- c(1:n)
  num <- seq(0, n, quotient)
  label <- c(1:K)
  
  for (i in 1:K){
    y[(num[i]+1):num[i+1]] <- label[i]
  }
  ##== Labeling End ==##
  
  # spatio temporal predictor
  n.t <- length(t) 
  
  # covariance matrix of gaussian process
  if (cov=="I") {
    Sigma <- error*diag(n.t)
  } else if(cov=="AR") {
    rho <- rho # 0.7
    idx.t <- c(1:n.t)
    Sigma <- error * outer(idx.t, idx.t, function(a, b) rho^abs(a-b))
  } else if(cov=="CS") {
    rho <- rho # 0.3
    tmp <- matrix(rho,n.t,n.t)
    diag(tmp) <- 1
    Sigma <- error * tmp
  } else warning("covariance structure was not specified")
  
  # empty x.list
  x.list <- as.list(1:n)
  
  # Divide index
  idx1 <- which(y==1)
  idx2 <- which(y==2)
  idx3 <- which(y==3)
  
  # Create mean vector
  mu.t1 <- 2*t
  mu.t2 <- -2*t
  mu.t3 <- rep(0.5,n.t)
  
  # x.list
  x.list <- as.list(1:n)
  for (i in 1:n) {
    if(i %in% idx1) x.list[[i]] <- rmvnorm(1, mu.t1, Sigma)
    else if(i %in% idx2) x.list[[i]] <- rmvnorm(1, mu.t2, Sigma)
    else if(i %in% idx3) x.list[[i]] <- rmvnorm(1, mu.t3, Sigma)
  }
  
  # Calculate the True p
  true.p <- rep(0,n)
  for(i in 1:n){
    a1 <- dmvnorm(x=x.list[[i]], mean = mu.t1,log=TRUE)
    a2 <- dmvnorm(x=x.list[[i]], mean = mu.t2,log=TRUE)
    a3 <- dmvnorm(x=x.list[[i]], mean = mu.t3,log=TRUE)
    total <- exp(a1)+exp(a2)+exp(a3)
    if (i %in% idx1) true.p[i] <- exp(a1)/total
    else if (i %in% idx2) true.p[i] <- exp(a2)/total
    else if (i %in% idx3) true.p[i] <- exp(a3)/total
    else warning(paste0("No class was assigned for obs ",i))
  }
  
  obj <- list(x = x.list, y = y, t = t, true.p = true.p)
  return(obj)
}
