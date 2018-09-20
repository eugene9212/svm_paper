pi.path.for.train<- 
  function(lambda, x.list, y, t, L, eps = 1.0e-5, Nmoves = 10 * n, ridge = 0)
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
    coef <- lapply(x.list, coef.fn, basis.mat = basis.mat) # list of 100(# of sample), ¾È¿¡´Â 10(# of basis ft)
    
    Phi <- inprod(basis.obj, basis.obj) # inprod(inner product of functional data objects)
    
    K <- matrix(0, n, n)
    
    # construct K
    for (i in 1:n) { 
      for (j in 1:i) {
        c1 <- as.matrix(coef[[i]])
        c2 <- as.matrix(coef[[j]])
        K[i,j] <- K[j,i] <- c1 %*% Phi %*% t(c2)# %*% rep(1, L)) %*% weight 
      }
    }
    
    # calculate pi
    p <- 0.4
    
    part1 <- fsvm.sub.pi.path(lambda,x.list, y, pi0 = p, K, eps, Nmoves, ridge = ridge)
    part2 <- fsvm.sub.pi.path(lambda,x.list,-y, pi0 = 1-p, K, eps, Nmoves, ridge = ridge)
    
    pi1 <- part1$Pi;          pi2 <- 1- part2$Pi
    n1 <- length(pi1);        n2 <- length(pi2)
    alpha01 <- part1$alpha0;  alpha02 <- -part2$alpha0
    alpha1 <- part1$alpha;    alpha2 <- part2$alpha
    
    alpha1 <- matrix(alpha1, length(alpha1)/length(pi1), length(pi1))
    alpha2 <- matrix(alpha2, length(alpha2)/length(pi2), length(pi2))
    
    result1 <- matrix(rbind(pi1, alpha01, alpha1), n+2, n1)
    result2 <- matrix(rbind(pi2, alpha02, alpha2), n+2, n2)
    result2 <- matrix(result2[,order(result2[1,])], n+2, n2)
    
    result <- cbind(result2[,-n2], result1[,-1])
    
    pi <- c(0, result[1,], 1)
    alpha0 <- c(result[2,1], result[2,], result[2,length(result[2,])])
    alpha <- cbind(rep(0, n), result[-(1:2),], rep(0, n))
    
    # calculate fx and gx
    alpha0 <- matrix(alpha0, n, length(alpha0) ,byrow = T)
    gx <- (K %*% (alpha * y))
    fx <- (alpha0 + gx)/lambda

    # calculate probability
    pi.star <- rep(0, new.n)
    
    for (i in 1:n) {
      minus <- which(sign(fx[i,])<0)
      if(min(minus) == 1 && (minus[2]-minus[1]) != 1){
        index1 <- minus[2]
      }else{
        index1 <- min(minus)
      }
      
      plus <- which(sign(fx[i,])>0)
      if(length(plus) == 0){
        pi.star[i] <- pi[index1] + 1.0e-8
      }else{
        index2 <- max(plus)
        
        pi.star[i] <- (pi[index1] + pi[index2])/2      
      }
    }
    
    # return
    obj <- list(lambda = lambda, pi = pi, prob = pi.star, alpha0  = alpha0, alpha = alpha, 
                fx = fx, gx = gx , L = L, t = t, k = K, basis.obj = basis.obj, coef = coef, 
                x = x.list, y = y)
  }
