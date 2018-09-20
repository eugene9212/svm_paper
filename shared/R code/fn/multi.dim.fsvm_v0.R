multi.dim.fsvm <- function(y, x.list, L = 10, lambda = 1, rho = 1, 
                           max.dist = sqrt(2) + 1.0e-4, weight = rep(1, n)) 
{
  n <- length(y)
  
  basis.obj <- coef <- as.list(1:n)
  for (i in 1:n)
  {
    xx <- x[[i]]
    rangeval <- quantile(t, c(0,1))
    basis.obj[[i]] <- create.bspline.basis(rangeval, L) 
    basis.mat <- eval.basis(t, basis.obj[[i]])
    
    cc <- matrix(0, n.grid, L)
    for (s in 1:n.grid) {
      cc[s,] <- lsfit(basis.mat, xx[s,], intercept = F)$coef  
    }
    coef[[i]] <- cc
  }
  
  # neighbor index
  neighbor.dist <- neighbor.index <- as.list(1:n.grid)
  s.dist <- as.matrix(dist(spatial.index, diag = T, upper = T))
  for (s in 1:n.grid) {
    neighbor.index[[s]] <- which(s.dist[s,] < max.dist)
    neighbor.dist[[s]] <- s.dist[s,neighbor.index[[s]]]
  }
  
  K <- matrix(0, n, n)
  for (i in 1:n)
  { 
    for (j in 1:n) 
    {
      Phi <- inprod(basis.obj[[i]], basis.obj[[j]])
      c1 <- coef[[i]]
      c2 <- coef[[j]]
      for (s in 1:n.grid) 
      {
        temp <- 0
        s.dist <- neighbor.dist[[s]]
        count <- 0
        for (u in neighbor.index[[s]]) {
          count <- count + 1
          dd <- s.dist[count]
          omega <- ifelse( dd < 1.0e-10, 1, rho^dd)
          temp <- temp + omega * t(c1[s,]) %*% Phi %*% c2[u,]  
        }
        K[i,j] <- K[i,j] + temp
      }
    }
  }
  
  # solve SVM
  Kscript <- (K * outer(y,y))
  cvec <- -rep(1,n)
  alpha <- wsvmQP(Kscript/lambda, cvec, weight, y)$alpha
  
  Elbow <- which(alpha > 1.0e-3 & alpha < 1 - 1.0e-3)
  gx <- (K %*% (alpha * y))
  if (length(Elbow) == 0) {
    alpha0 <- - mean(gx)
  } else {
    alpha0 <- mean(lambda * y[Elbow] - gx[Elbow])  
  }
  
  
  fx <- (alpha0 + gx)/lambda
  
  obj <- list(alpha = alpha, alpha0 = alpha0, fx = fx, 
              y = y, x = x, t = t, spatial.index = spatial.index,
              rho = rho, lambda = lambda, L = L, 
              neighbor.index = neighbor.index, 
              neighbor.dist = neighbor.dist,
              K = K, basis.obj = basis.obj, coef = coef)
}
