predict.multi.dim.fsvm <- function(obj, new.x) {
  new.n <- length(new.x)
  
  lambda <- obj$lambda
  alpha <- obj$alpha
  alpha0 <- obj$alpha0
  
  spatial.index <- obj$spatial.index
  n.grid <- nrow(spatial.index)
  t <- obj$t
  n.t <- length(t)
  
  y <- obj$y
  
  L <- obj$L
  coef <- obj$coef
  basis.obj <- obj$basis.obj
  
  neighbor.dist  <- obj$neighbor.dist
  neighbor.index <- obj$neighbor.index
  rho <- obj$rho
  
  
  new.basis.obj <- new.coef <- as.list(1:new.n)
  for (i in 1:new.n)
  {
    new.xx <- new.x[[i]]
    rangeval <- quantile(t, c(0,1))
    new.basis.obj[[i]] <- create.bspline.basis(rangeval, L) 
    new.basis.mat <- eval.basis(t, new.basis.obj[[i]])
    
    new.cc <- matrix(0, n.grid, L)
    for (s in 1:n.grid) {
      new.cc[s,] <- lsfit(new.basis.mat, new.xx[s,], intercept = F)$coef  
    }
    new.coef[[i]] <- new.cc
  }
  
  
  # compute K (new.n by n)
  new.K <- matrix(0, new.n, n)
  for (i in 1:new.n)
  { 
    for (j in 1:n) 
    {
      Phi <- inprod(new.basis.obj[[i]], basis.obj[[j]])
      c1 <- new.coef[[i]]
      c2 <- coef[[j]]
      for (s in 1:n.grid) 
      {
        temp <- 0
        s.dist <- neighbor.dist[[s]]
        count <- 0
        for (u in neighbor.index[[s]]) {
          count <- count + 1
          dd <- s.dist[count]
          omega <- ifelse(dd < 1.0e-10, 1, rho^dd)
          temp <- temp + omega * t(c1[s,]) %*% Phi %*% c2[u,]  
        }
        new.K[i,j] <- new.K[i,j] + temp
      }
    }
  }
  
  
  new.gx <- (new.K %*% (alpha * y))
  new.fx <- (alpha0 + new.gx)/lambda
  return(new.fx)
}
