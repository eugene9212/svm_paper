model2 <- function(n, n.basis, period, t = seq(0, 1, by = 0.05), n.r = 4, n.c = 4, seed = 1)
{
  set.seed(seed)
  # binary response
  y <- rep(1, n)
  y[1:(n/2)] <- -1
  
  # spatio temporal predictor
  n.t <- length(t) 
  
  n.grid <- n.r * n.c # of total lattice
  
  # spatial index matrix (r, c)
  spatial.index <- matrix(0, n.grid, 2)
  spatial.index[,1] <- mapply(rep, 1:n.r, n.c)
  spatial.index[,2] <- 1:n.c 
  colnames(spatial.index) <- c("row", "col")
  
  # neighbor index
  neighbor.index <- as.list(1:n.grid)
  s.dist <- as.matrix(dist(spatial.index, diag = T, upper = T))
  max.dist <- sqrt(2) + 1.0e-8
  for (s in 1:n.grid) {
    neighbor.index[[s]] <- which(s.dist[s,] < max.dist)
  }
  
  # create data
  basisobj <- create.fourier.basis(range(t), n.basis, period) 
  # plot(basisobj)
  
  D <- c(runif(n*n.basis, 0, 2)) 
  coefmat <- matrix(D, k, n)
  
  sbasis <- fd(coefmat, basisobj)
  tbasis <- fd(coefmat, basisobj)
  
  # correlation matrix
  Sigma <- 0.1 * outer(seq(0, n.basis, by = 1), seq(0, n.basis, by = 1), function(a, b) exp(-abs(a-b)))
  
  corrmat <- Sigma
  
  sbasis <- basisobj 
  tbasis <- basisobj 
  
  corrfd <- bifd(corrmat, sbasis, tbasis)
  
  
  # plot(tempfd)

  u.list <- as.list(1:n)
  for (i in 1:n) {
    mu.t <- z.t <- matrix(0, n.grid, n.t)
    # underlying covariate process, Z(t)
    for (s in 1:n.grid) {
      # mean vector of gaussian process
      mu.t[s,] <- sin(t)
      z.t[s,] <- rmvnorm(1, mu.t[s,], Sigma)
    }
    
    # spatially correlated process, U(t) = sum_neighbor Z(t)
    u.t <- matrix(0, n.grid, n.t)
    for (s in 1:n.grid) {
      u.t[s,] <- apply(z.t[neighbor.index[[s]],,drop = F], 2, mean)
    }
    u.list[[i]] <- u.t
  }

  x.list <- u.list
  for (i in 1:n) {
    if (y[i] == 1) {
      x.list[[i]] <-  u.list[[i]] + beta
    } else x.list[[i]] <- u.list[[i]]
  }
  obj <- list(x = x.list, y = y, t = t, spatial.index = spatial.index)
  return(obj)
}