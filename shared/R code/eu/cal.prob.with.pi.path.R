# This code calculate the probability with trained pi.path

cal.prob.with.pi.path <- function(obj, new.x) {
  
  new.n <- length(new.x)
  
  y <- obj$y ; n <- length(y)
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
  
  # # 0725ㅔ 추가해본거. fx가 infinity나오는거 같아서
  # new.fx1 <- matrix(0, dim(new.fx)[1], dim(new.fx)[2])
  # for (i in 1:dim(new.fx)[1])
  # { 
  #   for (j in 1:dim(new.fx)[2]) 
  #   {
  #     # if (new.fx[i,j] == Inf) {
  #     #   new.fx1[i,j] <- 10000
  #     # }else if(new.fx[i,j] == -Inf){
  #     #   new.fx1[i,j] <- -10000
  #     # }else 
  #       tryCatch(length(new.fx[i,j]) == 0, error = function(e) print("new.fx has length 0"))
  #     # }else{
  #     #   new.fx1[i,j] <- new.fx[i,j]
  #   }
  # }
  
  
  # calculate the probability
  pi.star <- rep(0, new.n)
  
  # pi.star
  for (i in 1:new.n) {
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
  
  # return(new.fx)
  obj <- list(prob = pi.star, new.fx = new.fx, new.K = new.K)
}
