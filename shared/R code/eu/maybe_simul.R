
########################################
## Simulation with train and test set ## (sin and cos)
########################################

rm(list = ls())

#### load packages
library(mvtnorm)
library(fda)

#### load R codes
setwd('C:/Users/eugene/Desktop/SVM_R/shared/R code/')
source('eu/sc.R')
source('eu/pi.path.for.train.R')
source('eu/cal.prob.with.pi.pathv0.R')

# source('fn/fsvm.pi.path.R')
source('fn/fsvm.sub.pi.path.R')
dyn.load("KernSurf/temp/wsvmqp.dll")

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('KernSurf/R')

####================================= <Simulation1> =====================================
# set up
# n.sim <- 100
n <- 10
t <- seq(0, 2*pi, by = 0.05)
lambda <- 0.9 # fsvm & pi.path
L <- 10

# storage
# result <- matrix(0, n.sim, 1)

# for (iter in 1:n.sim) {
#   tic <- Sys.time()
  iter <- 1
  
  # Data generation
  set.seed(iter)
  data <- sc(n, t, seed = iter)
  
  # train and test set
  id1 <- sample(n/2, n/4)
  id2 <- sample((n/2+1):n, n/4)
  id <- c(id1, id2)
  
  train.x <- data$x[id]
  train.y <- data$y[id]
  test.x <- data$x[-id]
  test.y <- data$y[-id]
  
  # plot train.x
  train.y
  plot(train.x[[1]]) # for y=-1
  plot(train.x[[3]]) # for y=+1
  
  obj_pi <- pi.path.for.train(lambda, train.x, train.y, t, L)
  
  # PLOT1
  pi <- obj_pi$pi
  alpha <- obj_pi$alpha
  plot(0,0, type = "l", xlim = c(0,1), ylim = c(0,1),
       xlab = "pi", ylab = "alpha")
  axis(1, round(pi,2), las=2)
  abline(v = pi, lty = 2, col = "gray")
  for (ii in 1:length(train.y)) lines(pi, alpha[ii,], col = 2)

  # PLOT2 (fx, also piece-wise linear)
  alpha <- obj_pi$alpha
  y <- obj_pi$y
  K <- obj_pi$K
  alpha0 <- matrix(obj_pi$alpha0, dim(gx)[1], length(obj_pi$alpha0) ,byrow = T)
  alpha0 <- matrix(alpha0, dim(gx)[1], dim(gx)[2], byrow = TRUE) 
  fx <- (alpha0 + gx)/lambda
  
  plot(0,0, type = "l", xlim = c(0,1), ylim = c(min(fx),max(fx)),
       xlab = "pi", ylab = "fx")
  axis(1, round(pi,2), las=2)
  abline(v = pi, lty = 2, col = "gray")
  for (ii in 1:length(y)) lines(pi, fx[ii,], col = 2)
  
  # calculate probability
  obj <- cal.prob.with.pi.path(obj_pi, test.x)
  
  # PLOT2 (fx, also piece-wise linear)
  fx <- obj$new.fx
  plot(0,0, type = "l", xlim = c(0,1), ylim = c(min(fx),max(fx)),
       xlab = "pi", ylab = "fx")
  axis(1, round(pi,2), las=2)
  abline(v = pi, lty = 2, col = "gray")
  for (ii in 1:length(test.y)) lines(pi, fx[ii,], col = 2)
  
  
  # calculate the probability
  pi.star <- rep(0, new.n)

  # pi.star
  for (i in 1:new.n) {
    minus <- which(sign(new.fx[i,])<0)
    if(min(minus) == 1 && (minus[2]-minus[1]) != 1){
      index1 <- minus[2]
    }else{
      index1 <- min(minus)
    }
    
    plus <- which(sign(new.fx[i,])>0)
    if(length(plus) == 0){
      pi.star[i] <- pi[index1] + 1.0e-8
    }else{
      index2 <- max(plus)
      
      pi.star[i] <- (pi[index1] + pi[index2])/2      
    }
  }
  # Boxplot of pi.star 
  # png(filename = paste0(iter,".png"))
  boxplot(pi.star[test.y == 1], pi.star[test.y != 1], xlab=iter)
  
  # Criteria
  CRE <- -1/length(test.y)*(sum(1/2*(1+test.y)*log(hat.pi) + 1/2*(1-test.y)*log(1-hat.pi)))
  boxplot(hat.pi[test.y == 1], hat.pi[test.y != 1], xlab=iter)
  
  # results  
  result[iter,] <- c(CRE)
  
  toc <- Sys.time()
  print(toc - tic)
}

result
# save results
write(result, "eu/result/result.txt")

####================================= <Simulation2> =====================================
# set up
n.sim <- 100
n <- 10
t <- seq(0, 1, by = 0.05)
lambda <- 0.5 # fsvm & pi.path
L <- 10

# storage
result <- matrix(0, n.sim, 1)

for (iter in 1:n.sim) {
  tic <- Sys.time()
  iter <- 1
  
  # Data generation
  set.seed(iter)
  data <- gp.1dim.sc(n, t, seed = iter)
  
  # train and test set
  id1 <- sample(n/2, n/4)
  id2 <- sample((n/2+1):n, n/4)
  id <- c(id1, id2)
  
  train.x <- data$x[id]
  train.y <- data$y[id]
  test.x <- data$x[-id]
  test.y <- data$y[-id]
  
  # plot train.x
  train.y
  plot(train.x[[1]]) # for y=-1
  plot(train.x[[3]]) # for y=+1
  
  print(iter)
  
  obj_pi <- pi.path.for.train(lambda, train.x, train.y, t, L)

  # PLOT2 (fx, also piece-wise linear)
  gx <- (obj_pi$K %*% (obj_pi$alpha * obj_pi$y))
  alpha0 <- matrix(obj_pi$alpha0, dim(gx)[1], length(obj_pi$alpha0) ,byrow = T)
  fx <- (alpha0 + gx)/obj_pi$lambda
  plot(0,0, type = "l", xlim = c(0,1), ylim = c(-10,10),
       xlab = "pi", ylab = "fx")
  abline(v = pi, lty = 2, col = "gray")
  for (ii in 1:length(train.y)) lines(obj_pi$pi, fx[ii,], col = 2)
  
  # calculate probability
  obj <- cal.prob.with.pi.path(obj_pi, test.x)
  # PLOT2 (fx, also piece-wise linear)
  fx <- obj$new.fx
  plot(0,0, type = "l", xlim = c(0,1), ylim = c(min(fx),max(fx)),
       xlab = "pi", ylab = "fx")
  axis(1, round(pi,2), las=2)
  abline(v = pi, lty = 2, col = "gray")
  for (ii in 1:length(test.y)) lines(pi, fx[ii,], col = 2)
  
  hat.pi <- obj$prob
  
  # Criteria  
        
  CRE <- -1/length(test.y)*(sum(1/2*(1+test.y)*log(hat.pi) + 1/2*(1-test.y)*log(1-hat.pi)))
  boxplot(hat.pi[test.y == 1], hat.pi[test.y != 1], xlab=iter)
  
  # results  
  result[iter,] <- c(CRE)
  
  toc <- Sys.time()
  print(toc - tic)
}

result
# save results
write(result, "eu/result/result.txt")


####================================= <Simulation3> =====================================

source('eu/gp.1dim.R')

# set up
n.sim <- 100
n <- 10
beta <- 1
t <- seq(0, 1, by = 0.05)
lambda <- 0.5 # fsvm & pi.path
L <-10

# storage
result <- matrix(0, n.sim, 1)

for (iter in 1:n.sim) {
  tic <- Sys.time()
  iter <- 1
  
  # Data generation
  set.seed(iter)
  data <- gp.1dim(n, t, seed = iter)
  
  # train and test set
  id1 <- sample(n/2, n/4)
  id2 <- sample((n/2+1):n, n/4)
  id <- c(id1, id2)
  
  train.x <- data$x[id]
  train.y <- data$y[id]
  test.x <- data$x[-id]
  test.y <- data$y[-id]
  
  print(iter)
  
  obj_pi <- pi.path.for.train(lambda, train.x, train.y, t, L)
  # PLOT1
  pi <- obj_pi$pi
  alpha <- obj_pi$alpha
  plot(0,0, type = "l", xlim = c(0,1), ylim = c(0,1),
       xlab = "pi", ylab = "alpha")
  axis(1, round(pi,2), las=2)
  abline(v = pi, lty = 2, col = "gray")
  for (ii in 1:length(train.y)) lines(pi, alpha[ii,], col = 2)

  # PLOT2 (fx, also piece-wise linear)
  gx <- (obj_pi$K %*% (obj_pi$alpha * obj_pi$y))
  alpha0 <- matrix(obj_pi$alpha0, dim(gx)[1], length(obj_pi$alpha0) ,byrow = T)
  fx <- (alpha0 + gx)/obj_pi$lambda
  plot(0,0, type = "l", xlim = c(0,1), ylim = c(min(fx),max(fx)),
       xlab = "pi", ylab = "fx")
  axis(1, round(pi,2), las=2)
  abline(v = pi, lty = 2, col = "gray")
  for (ii in 1:length(train.y)) lines(pi, fx[ii,], col = 2)
  
  # calculate probability
  obj <- cal.prob.with.pi.path(obj_pi, test.x)
  # PLOT2 (fx, also piece-wise linear)
  fx <- obj$new.fx
  plot(0,0, type = "l", xlim = c(0,1), ylim = c(min(fx),max(fx)),
       xlab = "pi", ylab = "fx")
  axis(1, round(pi,2), las=2)
  abline(v = pi, lty = 2, col = "gray")
  for (ii in 1:length(test.y)) lines(pi, fx[ii,], col = 2)
  
  hat.pi <- obj$prob
  
  # Criteria
  CRE <- -1/length(test.y)*(sum(1/2*(1+test.y)*log(hat.pi) + 1/2*(1-test.y)*log(1-hat.pi)))
  boxplot(hat.pi[test.y == 1], hat.pi[test.y != 1], xlab=iter)
  
  # results  
  result[iter,] <- c(CRE)
  
  toc <- Sys.time()
  print(toc - tic)
}

result
# save results
write(result, "eu/result/result.txt")


