"predict.fsvm.pi.path" <- function(obj, new.x) {
  
  test.n <- length(new.x)
  
  y <- obj$y ; train.n <- length(y)
  k <- obj$k
  lambda <- obj$lambda
  alpha <- obj$alpha
  alpha0 <- obj$alpha0
  pi <- obj$pi

  t <- obj$t
  n.t <- length(t)
  
  L <- obj$L
  coef <- obj$coef
  basis.obj <- obj$basis.obj
  
  basis.mat <- eval.basis(t, basis.obj)
  
  coef.fn <- function(xx, basis.mat) {
    d.mat <- t(basis.mat) %*% basis.mat
    chol.obj <- chol(d.mat)
    coef <- xx %*% basis.mat %*% chol2inv(chol.obj) # ls
    return(coef)
  }
  new.coef <- lapply(test.x, coef.fn, basis.mat = basis.mat) # Easily parellelizable!
  
  # kernel matrix for new.x
  Phi <- inprod(basis.obj, basis.obj)
  
  new.K <- matrix(0, test.n, train.n)
  for (i in 1:test.n)
  { 
    for (j in 1:train.n) 
    {
      c1 <- as.matrix(new.coef[[i]])
      c2 <- as.matrix(coef[[j]])
      new.K[i,j] <- c1 %*% Phi %*% t(c2) # %*% rep(1, L)  * weight # 여기 필요없는거 맞나??
    }
  }
  
  alpha0 <- matrix(alpha0, test.n, length(alpha0) ,byrow = T)
  new.gx <- (new.K %*% (alpha * y))
  new.fx <- (alpha0 + new.gx)/lambda
  
  # calculate the probability
  pi.star <- rep(0, test.n)
  
  # pi.star
  for (i in 1:test.n) {
    minus <- which(sign(new.fx[i,])<0)
    if(min(minus) == 1 && (minus[2]-minus[1]) != 1){
      index1 <- minus[2]
    }else{
      index1 <- min(minus)
    }
    
    plus <- which(sign(new.fx[i,])>0)
    if(length(plus) == 0){
      pi.star[i] <- pi[index1] + 1.0e-8
    }else{
      index2 <- max(plus)
      
      pi.star[i] <- (pi[index1] + pi[index2])/2      
    }
  }
  
  obj <- list(prob = pi.star, new.fx = new.fx, new.K = new.K)
}
