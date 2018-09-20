"fsvm.pi.path" <- 
  function(lambda, x, y, K,  
           eps = 1.0e-5, Nmoves = 5 * n, ridge = 0)
  {
    n <- length(y)
    p <- 0.5
    
    part1 <- fsvm.sub.pi.path(lambda,x, y, pi0 = p, K, eps, Nmoves, ridge = ridge)
    part2 <- fsvm.sub.pi.path(lambda,x,-y, pi0 = 1-p, K, eps, Nmoves, ridge = ridge)
    
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
    
    obj <- list(lambda = lambda, pi = pi, alpha0  = alpha0, alpha = alpha, 
                x = x, y = y)
    obj
  }
