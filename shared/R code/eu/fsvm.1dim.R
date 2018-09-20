# what is neightbor.index? 
# 지워도 돼는건지? 일차원에선 신경안써도 돼는건지?
# example
# rm(list = ls())
# n <- 100 ; beta <- 1
# a <- gp.1dim(n, beta, t = seq(0, 1, by = 0.05), seed = 10) ; x.list <- x <- a$x ;y <- a$y
# t = seq(0, 1, by = 0.05)
# solve one dimensional fsvm

fsvm.1dim <- function(y, x.list, t, L = 10, lambda = 1, rho = 1, weight = rep(1, n)) 
{
  n <- length(y)
  
  # basis matrix for a common time grid
  rangeval <- quantile(t, c(0,1))
  basis.obj <- create.bspline.basis(rangeval, L) # L : number of basis
  basis.mat <- eval.basis(t, basis.obj) # 21(# of t)*10(# of basis ft)
  
  # compute coefficient (least square fit)
  coef.fn <- function(xx, basis.mat) {
    d.mat <- t(basis.mat) %*% basis.mat
    chol.obj <- chol(d.mat)
    coef <- xx %*% basis.mat %*% chol2inv(chol.obj) # ls
    return(coef)
  }
  coef <- lapply(x.list, coef.fn, basis.mat = basis.mat) # list of 100(# of sample), 안에는 10(# of basis ft)

  Phi <- inprod(basis.obj, basis.obj) # 10*10 / inprod(inner product of functional data objects)
  weight <- rho
  
  K <- matrix(0, n, n)
  for (i in 1:n) { 
    for (j in 1:i) {
      c1 <- coef[[i]]
      c2 <- coef[[j]]
      K[i,j] <- K[j,i] <- t(((c1[1,] %*% Phi) * c2[1,,drop = F]) %*% rep(1, L)) %*% weight 
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
              y = y, x = x.list, t = t,
              rho = rho, lambda = lambda, L = L,
              K = K, basis.obj = basis.obj, coef = coef)
}

# result <- fsvm.1dim(y, x.list, L = 10, lambda = 1, rho = 1, weight = rep(1, n))
