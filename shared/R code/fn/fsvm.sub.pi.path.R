"fsvm.sub.pi.path" <- 
  function(lambda, x, y, 
           pi0 = 0.5, K, 
           eps = 1e-5, Nmoves = 100 * n, ridge = 0)
  {
    n <- length(y)
    index.plus <- which(y == 1)
    index.minus <- which(y == -1)
    nplus <- length(index.plus)
    nminus <- length(index.minus)
    Right <- Elbow <- NULL
    Left <- seq(n)
    Kscript <- K * outer(y, y)
    Kstar <- matrix(0, 1, 1)
    
    alpha <- matrix(0, n, Nmoves)      # alpha matrix (Nmoves*n)
    alpha0 <- double(Nmoves)           # alpha0 vector
    Elbow.list <- as.list(seq(Nmoves)) # Elbow storage list
    Left.list <- as.list(seq(Nmoves))  # Elbow storage list
    Pi <- double(Nmoves)               # Pi vector
    
    #
    # Initialization for alpha/alpha0/Sets (given lambda, and pi0)
    #
    init <- pi.initialization(pi0, x, y, K, lambda) # subroutine
    
    # from 'init' object
    Elbow.list[[1]]<-Elbow<-init$Elbow # Elbow 
    Left.list[[1]]<-Left <- init$Left  # Left  
    Right <- init$Right                # Right 
    alpha0[1] <- init$alpha0           # alpha0
    alpha[,1] <- init$alpha            # alpha 
    Pi[1] <- pi0                       # pi    
    
    # Updata Kstar from the initialized values (init)
    first <- c(0, y[Elbow])                               # first row
    rest <- matrix(cbind(y[Elbow], Kscript[Elbow, Elbow]),# rest rows
                   length(Elbow), length(Elbow) + 1)      # (Elbow * Elbow+1)
    Kstar <- rbind(first, rest)                           # combine
    fl <- (K %*% (alpha[,1] * y) + alpha0[1])/lambda      # fl
    k <- 1
    
    for (k in 1:Nmoves){
      if (length(Elbow) == 0) {
        nlplus <- length(intersect(Left, index.plus))
        nl <- length(Left)
        balanced <- abs(nlplus/nl - Pi[k]) < eps
        if (balanced){
          empty <- wsvm.empty.pi.path(pi = Pi[k], a = alpha[,k], a0 = alpha[k],
                                      lambda, x, y, K, Left, Right)
          Pi[k+1] <- empty$pi
          alpha[,k+1] <- empty$alpha
          alpha0[k+1] <- empty$alpha0
          Elbow <- empty$Elbow
          Left <- empty$Left
          Right <- empty$Right
          
          first <- c(0,y[Elbow])
          rest <- matrix(cbind(y[Elbow], Kscript[Elbow, Elbow]), 1, 2)
          Kstar <- rbind(first, rest)
          fl <- (K %*% (alpha[,k+1]*y) + alpha0[k+1])/lambda
        } else stop("Empty elbow situation with unbalanced left set")
      }   
      else {
        ne <- length(Elbow)# n_elbow
        nl <- length(Left) # n_left
        Kle <- matrix(K[Left,Elbow], nl, ne) # compute K(nl*ne)
        kl <- y[Elbow] * apply(Kle , 2, sum) # compute kl
        bstar <- PiSolveKstar(Kstar, vec = kl, con = nl, ridge = ridge) # solve linear equation
        b0 <- bstar[1] 
        b <- bstar[-1] 
        gl <- (matrix(K[,Elbow],n,ne) %*% (y[Elbow]*b)-
                 matrix(K[,Left],n,nl) %*% rep(1,nl) + b0) / lambda
        temp <- -alpha[Elbow, k] + Pi[k]*b               
        pi.left <- (y[Elbow]==-1)*temp/(b-1) +           
          (y[Elbow]==1)*(1+temp)/(1+b)          # ** tricky **(different for y_i)
        pi.right <- temp/b                               # exit to Right (alpha = 0)
        pi01 <- c(pi.right, pi.left)                     # combine
        pi.exit <- min(pi01[pi01 > (Pi[k] + eps)],2)     # find smallest point
        pi.exit[is.nan(pi.exit)+is.na(pi.exit) != 0] <- 2# exclue NaN/NA
        
        pii <- (y - fl)/gl + Pi[k]   
        pii[is.nan(pii)] <- -2
        pi.entry <- min(pii[pii > (Pi[k] + eps)], 2)
        pi.entry[is.nan(pi.entry)+is.na(pi.entry) != 0] <- 2
        
        new.pi <- min(pi.entry, pi.exit)  # smaller(exit, enter)
        Pi[k + 1] <- new.pi               # update pi
        w <- wvec(Pi[k+1], y = y)         # update weight vector
        
        alpha[,k+1] <- alpha[,k] 
        alpha[Elbow,k+1] <- alpha[Elbow,k] - (Pi[k]-Pi[k+1])*b # update alpha in Elbow
        alpha[Left,k+1] <- w[Left]                             # update alpha in Left
        alpha0[k+1] <- alpha0[k] - (Pi[k]-Pi[k+1])*b0          # update alpha0
        fl <- fl + (Pi[k+1] - Pi[k]) * gl                      # update fl
        #if (abs(sum(fl) + n) < eps) {
        #  immobile = prod(abs(fl + 1) < eps)
        #  if (immobile == 1) break }
        
        # a) If "Exit from Elbow"
        if (pi.exit < pi.entry){ 
          i1 <- Elbow[pi.left == new.pi]
          i2 <- Elbow[pi.right == new.pi]
          i <- union(i1, i2)
          Elbow <- setdiff(Elbow, i)
          Left <- union(Left, i1)
          Right <- union(Right, i2)} 
        
        # b) If "Enter to Elbow"
        else {     
          i1 <- Left[pii[Left] == new.pi]
          i2 <- Right[pii[Right] == new.pi]
          Elbow <- c(Elbow, i1, i2)
          Left <- setdiff(Left, i1)
          Right <- setdiff(Right,i2)}
        
        first <- c(0, y[Elbow])
        rest <- matrix(cbind(y[Elbow], Kscript[Elbow, Elbow]), length(Elbow), length(Elbow) + 1)
        Kstar <- rbind(first, rest)}
      
      if (Pi[k+1] > 1-eps) break # terminate updating
      k <- k + 1                 # update k
      Elbow.list[[k]] <- Elbow   # save Elbow to the list
      Left.list[[k]] <- Left   # save Elbow to the list
    } # End of the update
    
    obj <- list(Pi = Pi[seq(k)], alpha = alpha[, seq(k)], alpha0 = alpha0[seq(k)], pi = Pi[seq(k)],
                Elbow = Elbow.list[seq(k)], Left = Left.list[seq(k)], pi0 = pi0,
                x = x, y = y, lambda = lambda)
    obj
  }
