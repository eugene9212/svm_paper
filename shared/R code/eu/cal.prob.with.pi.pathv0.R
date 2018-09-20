# This code calculate the probability with trained pi.path

cal.prob.with.pi.pathv0 <- function(obj, new.x) {
  # obj is from pi.path.for.train
  # obj <- list(lambda = lambda, pi = pi, alpha0  = alpha0, alpha = alpha, 
  #             L = L, t = t, K = K, basis.obj = basis.obj, coef = coef, 
  #             x = x.list, y = y)
  
  new.n <- length(new.x)
  
  y <- obj$y ; n <- length(y)
  lambda <- obj$lambda
  alpha <- obj$alpha
  alpha0 <- obj$alpha0
  pi <- obj$pi

  t <- obj$t
  n.t <- length(t)

  L <- obj$L
  coef <- obj$coef
  basis.obj <- obj$basis.obj
  
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
  
  new.K <- matrix(0, new.n, n)
  for (i in 1:new.n)
  { 
    for (j in 1:n) 
    {
      c1 <- as.matrix(new.coef[[i]])
      c2 <- as.matrix(coef[[j]])
      new.K[i,j] <- c1 %*% Phi %*% t(c2) # %*% rep(1, L)  * weight # 여기 필요없는거 맞나??
    }
  }
  alpha0 <- matrix(alpha0, new.n, length(alpha0) ,byrow = T)
  new.gx <- (new.K %*% (alpha * y))
  new.fx <- (alpha0 + new.gx)/lambda

  # return(new.fx)
  obj <- list(new.gx = new.gx, new.fx = new.fx, new.K = new.K, alpha = alpha, alpha0 = alpha0)
}
