"predict.fsvm.prob" <- function(obj, new.x) {
  
  test.n <- length(new.x)
  
  y <- obj$y ; train.n <- length(y)
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
  new.coef <- lapply(new.x, coef.fn, basis.mat = basis.mat) # Easily parellelizable!
  
  # kernel matrix for new.x
  Phi <- inprod(basis.obj, basis.obj)
  
  coef.m <- matrix(unlist(coef), nrow = length(coef), byrow=T)
  coef.m.new <- matrix(unlist(new.coef), nrow = length(new.coef), byrow=T)
  new.K <- coef.m.new %*% Phi %*% t(coef.m)

  alpha0 <- matrix(alpha0, test.n, length(alpha0), byrow = T)
  new.gx <- (new.K %*% (alpha * y))
  new.fx <- (alpha0 + new.gx)/lambda
  
  # calculate the probability
  pi.star <- rep(0, test.n)
  
  # revised(0912
  for (i in 1:test.n) {
    minus <- which(sign(new.fx[i,])<0)
    if (length(minus) == 0){
      pi.star[i] <- 1
      next
    # } else if (min(minus) == 1 && (minus[2]-minus[1]) != 1){
      # index1 <- minus[2]
    } else {
      index1 <- min(minus)
    }
    
    plus <- which(sign(new.fx[i,])>0)
    
    if(length(plus) == 0){
      pi.star[i] <- 0 #pi[index1] + 1.0e-8
    }else{
      index2 <- max(plus)
      
      pi.star[i] <- (pi[index1] + pi[index2])/2      
    }
  }
  
  obj <- list(prob = pi.star, new.fx = new.fx, new.K = new.K)
}
