fsvm.prob <- 
  function(x.list, y, t, L, eps = 1.0e-5, Nmoves = 100 * n, ridge = 0)
  {
    # x.list <- train.x ; y <- train.y
    n <- length(y)
    # eps = 1.0e-5; Nmoves = 100 * n; ridge = 0
    
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
    coef <- lapply(x.list, coef.fn, basis.mat = basis.mat) # list of 100(# of sample)
    
    Phi <- inprod(basis.obj, basis.obj) # inprod(inner product of functional data objects)
    
    # construct K
    coef.m <- matrix(unlist(coef), nrow = length(coef), byrow=T)
    K <- coef.m %*% Phi %*% t(coef.m)
    
    # Start lambda grid
    lambda <- c(0.1,0.5,0.9)
    z <- length(lambda)
    CRE <- matrix(0, z, 1)
    pi.list <- as.list(1:z)
    pi.star.list <- as.list(1:z)
    alpha0.list <- as.list(1:z)  
    alpha.list <- as.list(1:z)  
    fx.list <- as.list(1:z)
    gx.list <- as.list(1:z)
    
    # Loop
    for (ii in 1:z){
      lamb <- lambda[ii]
      
      # calculate pi
      p <- 0.4
      
      part1 <- fsvm.sub.pi.path(lamb,x.list, y, pi0 = p, K, eps, Nmoves, ridge = ridge)
      part2 <- fsvm.sub.pi.path(lamb,x.list,-y, pi0 = 1-p, K, eps, Nmoves, ridge = ridge)
      
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
      
      # (revised part : calculate the alpha0 when pi is 0 and 1)
      zero <- lamb
      one <- -lamb
      pi <- c(0, result[1,], 1)
      # alpha0 <- c(result[2,1], result[2,], result[2,length(result[2,])])
      alpha0 <- c(zero, result[2,], one)
      alpha <- cbind(rep(0, n), result[-(1:2),], rep(0, n))
      
      # calculate fx and gx
      alpha00 <- matrix(alpha0, n, length(alpha0), byrow = T)
      gx <- (K %*% (alpha * y))
      fx <- (alpha00 + gx)/lamb
      
      # calculate probability
      pi.star <- rep(0, n)
      
      # revised(0912)
      for (i in 1:n) {
        minus <- which(sign(fx[i,])<0)
        if (length(minus) == 0){
          pi.star[i] <- 1
          next
          # } else if (min(minus) == 1 && (minus[2]-minus[1]) != 1){
          # index1 <- minus[2]
        } else {
          index1 <- min(minus)
        }
        
        plus <- which(sign(fx[i,])>0)
        
        if(length(plus) == 0){
          pi.star[i] <- 0
        }else{
          index2 <- max(plus)
          
          pi.star[i] <- (pi[index1] + pi[index2])/2      
        }
      }
      
      # Criteria
      delta = 1.0e-8
      CRE[ii,] <- -1/length(y)*(sum(1/2*(1+y)*log(pi.star + delta) + 1/2*(1-y)*log(1-pi.star + delta)))
      pi.list[[ii]] <- pi
      pi.star.list[[ii]] <- pi.star
      alpha0.list[[ii]] <- alpha0
      alpha.list[[ii]] <- alpha
      fx.list[[ii]] <- fx
      gx.list[[ii]] <- gx
    }
    ind <- which(CRE == min(CRE))[1]
    opt.lambda <- lambda[ind]
    pi <- pi.list[[ind]]; pi.star <- pi.star.list[[ind]]; alpha0 <- alpha0.list[[ind]]
    alpha <- alpha.list[[ind]]; fx <- fx.list[[ind]]; gx <- gx.list[[ind]]
    
    # return
    obj <- list(lambda = opt.lambda, pi = pi, prob = pi.star, alpha0  = alpha0, alpha = alpha, 
                fx = fx, gx = gx , L = L, t = t, K = K, basis.obj = basis.obj, coef = coef, 
                x = x.list, y = y)
  }
