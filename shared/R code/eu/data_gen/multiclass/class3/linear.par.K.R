# create Gaussian Process
linear.par.K <- function(n, error, beta, K, cov, rho, t = seq(0, 1, by = 0.05), seed = 1)
{
  set.seed(seed)
  
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
    # create idx for each lable
    eval(parse(text=paste0(
      "idx",i," <- (num[",i,"]+1):num[",i,"+1]")))
  }
  ##== Labeling End ==##
  
  beta1 <- seq(0, (K-1), 1) * beta
  
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
  
  for (i in 1:K){
    
    # mean vector of gaussian process
    eval(parse(text=paste0(
      "mu.t",i,"<- t + beta1[",i,"]")))
    
    for(j in (num[i]+1):num[i+1]){
      
      # assign for each list
      eval(parse(text=paste0(
        "x.list[[",j,"]] <- rmvnorm(1, mu.t",i,", Sigma)")))
    }
  }
  
  # Calculate the True p
  true.p <- rep(0,n)
  
  for(i in 1:n){
    
    total <- 0
    
    for(j in 1:K){
      
    eval(parse(text=paste0(
      "a",j," <- dmvnorm(x=x.list[[",i,"]], mean = mu.t",j,",log=TRUE)")))
      
      for(jj in 1:K){
      eval(parse(text=paste0(
        "total <- c(total,exp(a",j,"))")))
      }
      
      total1 <- sum(total)
    }
     
    for(m in 1:K){
      eval(parse(text=paste0(
        "if (",i," %in% idx",m,") true.p[",i,"] <- exp(a",m,")/total1")))
    }
  }
  
  obj <- list(x = x.list, y = y, t = t, true.p = true.p)
  return(obj)
}
