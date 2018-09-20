multi.dim.fsvm <- function(y, x.list, spatial.index, L = 10, lambda = 1, rho = 1, 
                           max.dist = sqrt(2) + 1.0e-4, weight = rep(1, n)) 
{
  n <- length(y)
  n.grid <- nrow(spatial.index)
  
  # basis matrix for a common time grid
  rangeval <- quantile(t, c(0,1))
  basis.obj <- create.bspline.basis(rangeval, L) 
  basis.mat <- eval.basis(t, basis.obj)
  
  # compute coefficient (least square fit)
  coef.fn <- function(xx, basis.mat) {
    d.mat <- t(basis.mat) %*% basis.mat
    chol.obj <- chol(d.mat)
    coef <- xx %*% basis.mat %*% chol2inv(chol.obj) # ls
    return(coef)
  }
  x.list <- x # list of 20
  coef <- lapply(x.list, coef.fn, basis.mat = basis.mat) # Easily parellelizable!
  
  # neighbor index
  neighbor.dist <- neighbor.index <- as.list(1:n.grid)
  s.dist <- as.matrix(dist(spatial.index, diag = T, upper = T)) # this can be further improved via GPU computing
  
  
  for (s in 1:n.grid) {
    neighbor.index[[s]] <- which(s.dist[s,] < max.dist)
    neighbor.dist[[s]] <- s.dist[s,neighbor.index[[s]]]
  }
  
  
  Phi <- inprod(basis.obj, basis.obj)
  n.neighbor <- sapply(neighbor.index, length)
  index1 <- unlist(mapply(rep, 1:n.grid, n.neighbor))
  index2 <- unlist(neighbor.index)
  weight <- rho^unlist(neighbor.dist)
  
  K <- matrix(0, n, n)
  for (i in 1:n) { 
    for (j in 1:i) {
      c1 <- coef[[i]]
      c2 <- coef[[j]]
      K[i,j] <- K[j,i] <- t(((c1[index1,] %*% Phi) * c2[index2,,drop = F]) %*% rep(1, L)) %*% weight 
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
              s.dist = s.dist,
              K = K, basis.obj = basis.obj, coef = coef)
}
