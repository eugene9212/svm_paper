predict.1dim.fsvm <- function(obj, new.x) {
  # new.x <- test.x
  new.n <- length(new.x)
  
  lambda <- obj$lambda
  alpha <- obj$alpha
  alpha0 <- obj$alpha0
  
  t <- obj$t ; n.t <- length(t)
  y <- obj$y
  L <- obj$L
  coef <- obj$coef
  basis.obj <- obj$basis.obj
  
  rho <- obj$rho
  
  rangeval <- quantile(t, c(0,1))
  new.basis.obj <- create.bspline.basis(rangeval, L) 
  new.basis.mat <- eval.basis(t, new.basis.obj)
  
  coef.fn <- function(xx, basis.mat) {
    d.mat <- t(basis.mat) %*% basis.mat
    chol.obj <- chol(d.mat)
    coef <- xx %*% basis.mat %*% chol2inv(chol.obj) # ls
    return(coef)
  }
  new.coef <- lapply(new.x, coef.fn, basis.mat = new.basis.mat) # Easily parellelizable!
  
  # kernel matrix for new.x
  Phi <- inprod(basis.obj, basis.obj)
  weight <- rho
  
  n <- length(y)
  new.K <- matrix(0, new.n, n)
  for (i in 1:new.n)
  { 
    for (j in 1:n) 
    {
      c1 <- as.vector(new.coef[[i]])
      c2 <- as.vector(coef[[j]])
      new.K[i,j] <- (c1 %*% Phi %*% c2) # %*% rep(1, L)  * weight # 여기 필요없는거 맞나??
      # c1 <- new.coef[[i]]
      # c2 <- coef[[j]]
      # new.K[i,j] <- t(((c1[1,] %*% Phi) * c2[1,,drop = F]) %*% rep(1, L)) %*% weight 
    }
  }
  
  new.gx <- (new.K %*% (alpha * y))
  new.fx <- (alpha0 + new.gx)/lambda
  
  # return(new.fx)
  obj <- list(new.fx = new.fx, new.K = new.K)
}
