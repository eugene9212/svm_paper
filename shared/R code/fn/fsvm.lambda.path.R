"fsvm.lambda.path" <-
  function(w, x, y, K,
           eps = 1e-10, Nmoves = 5 * n, lambda.min = 1e-04, ridge = 0)
  {
    require(svmpath)
    n <- length(y)   
    yvals <- table(y)
    index.plus <- which(y == 1)
    index.minus <- which(y == -1)
    nplus <- yvals[2]
    nminus <- yvals[1]
    Right <- Elbow <- NULL
    Left <- seq(n)

    Kscript <- K * outer(y, y)
    alpha <- matrix(w, n, Nmoves)
    alpha0 <- double(Nmoves)
    Elbow.list <- as.list(seq(Nmoves)) 
    lambda <- double(Nmoves) 
    Kstar <- matrix(0, 1, 1) 
    
    if (abs(sum(w[index.plus]) - sum(w[index.minus])) < 1e-12) { 
      init <- wBalanced.Initialization(w, K, y, nplus, nminus)
      Elbow <- init$Elbow
      Left <- setdiff(Left, Elbow)
    }   else { 
      init <- wUnbalanced.Initialization(w, K, y, nplus, nminus)
      Elbow <- init$Elbow
      Right <- init$Right
      Left <- init$Left
    }
    Elbow.list[[1]] <- Elbow
    lambda0 <- init$lambda
    alpha00 <- init$alpha00
    alpha0[1] <- init$alpha0
    alpha[, 1] <- init$alpha
    Kstar <- UpdateKstar(Kstar, Kscript[Elbow, Elbow], NULL, y[Elbow])
    lambda[1] <- lambda0 
    fl <- (K %*% (alpha[, 1] * y) + init$alpha0)/lambda0
    
    k <- 1
    
    while ((k < Nmoves) && (lambda[k] > lambda.min)) {
      if (length(Elbow) == 0) {
        if (sum((y*w)[Left]) > eps)
          stop("Unbalanced data in interior empty elbow situation")
        init <- wBalanced.Initialization(w[Left], K[Left, Left], y[Left], 
                                         sum(y[Left]>0), sum(y[Left]<0))
        lambda0 <- init$lambda
        alpha0[k + 1] <- init$alpha0
        Elbow <- Left[init$Elbow]
        Left <- setdiff(Left, Elbow)
        lambda[k + 1] <- lambda0
        alpha[, k + 1] <- alpha[, k]
        Kstar <- UpdateKstar(Kstar, Kscript[Elbow, Elbow], 
                             NULL, y[Elbow])
        fl <- (lambda[k]/lambda[k + 1]) * (fl + (alpha0[k + 1] - alpha0[k])/lambda[k])
      } else {
        bstar <- SolveKstar(Kstar, ridge = ridge)
        b0 <- bstar[1]
        b <- bstar[-1]
        gl <- K[, Elbow, drop = FALSE] %*% (y[Elbow] * b) + b0
        dl <- fl - gl
        immobile <- sum(abs(dl))/n < eps
        temp <- -alpha[Elbow, k] + lambda[k] * b 
        lambda.left <- (w[Elbow] + temp)/b
        lambda.left[abs(b) < eps] <- -1
        lambda.right <- (0+temp)/b
        lambda.right[abs(b) < eps] <- -1
        lambda01 <- c(lambda.right, lambda.left)
        lambda.exit <- max(lambda01[lambda01 < lambda[k] - eps])
        
        if (immobile & (lambda.exit < eps)) {break
        } else if (!immobile) {
          lambdai <- (lambda[k] * (dl))/(y - gl)
          lambdai[abs(y - gl) < eps] <- -Inf     
          lambda.entry <- max(lambdai[lambdai < lambda[k] - eps])
        } else lambda.entry <- -1
        
        lambda.max <- max(lambda.entry, lambda.exit)
        lambda[k + 1] <- lambda.max
        alpha[, k + 1] <- alpha[, k]
        alpha[Elbow, k + 1] <- alpha[Elbow, k] - (lambda[k] - lambda[k + 1]) * b 
        alpha0[k + 1] <- alpha0[k] - (lambda[k] - lambda[k + 1]) * b0 
        fl <- (lambda[k]/lambda[k + 1]) * (dl) + gl
        if (lambda.entry > lambda.exit) {
          i <- match(lambda.entry, lambdai, 0)[1]
          if (match(i, Left, FALSE)) {
            Left <- setdiff(Left, i)
          } else {
            Right <- setdiff(Right, i)
          }
          Kstar <- UpdateKstar(Kstar, Kscript[i, i], drop(Kscript[i, Elbow]), y[i])
          Elbow <- c(Elbow, i)
        } else {
          idrop <- Leaveright <- NULL
          i <- Elbow[abs(lambda.right - lambda.exit) < eps]
          if (length(i) > 0) {
            Leaveright <- rep(TRUE, length(i))
            idrop <- i
          }
          i <- Elbow[abs(lambda.left - lambda.exit) < eps]
          if (length(i) > 0) {
            Leaveright <- c(Leaveright, rep(FALSE, length(i)))
            idrop <- c(idrop, i)
          }
          for (j in seq(along = idrop)) {
            if (Leaveright[j]) {
            } else {
              Left <- c(Left, idrop[j])
            }
            mi <- match(idrop[j], Elbow)
            Kstar <- DowndateKstar(Kstar, mi)
            Elbow <- Elbow[-mi]
          }
        }
      }
      k <- k + 1
      Elbow.list[[k]] <- Elbow
    }
    pi <- w[1]*(y[1] == -1) + (1-w[1])*(y[1] == 1)
    obj <- list(pi = pi, lambda = lambda[seq(k)], 
                alpha0 = alpha0[seq(k)], alpha = alpha[, seq(k)], x = x, y = y)
    obj
  }
