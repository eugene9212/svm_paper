# create K
# Create K instead multi.dim.fsvm (start) #
create.K <- function(x, t, L){
  
  n <- length(y)
  
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
  coef <- lapply(x, coef.fn, basis.mat = basis.mat) # Easily parellelizable!
  
  # calculate phi
  Phi <- inprod(basis.obj, basis.obj)
  
  K <- matrix(0, n, n)
  for (i in 1:n) { 
    for (j in 1:i) {
      c1 <- as.matrix(coef[[i]])
      c2 <- as.matrix(coef[[j]])
      K[i,j] <- K[j,i] <- c1 %*% Phi %*% t(c2) 
    }
  }
  
  obj <- list(K = K)
}