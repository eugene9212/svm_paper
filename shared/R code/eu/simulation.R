rm(list = ls())
# load packages
library(mvtnorm)
library(fda)

# load R codes
setwd('C:/Users/eugene/Desktop/SVM_R/shared/R code/')
source('eu/fsvm.1dim.R')
source('eu/fsvm.1dim.fourier.R')
source('eu/gp.1dim.R')
source('fn/fsvm.pi.path.R')
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

####================================= <Simulation 1> =====================================
# set up
n.sim <- 100
n <- 100
beta <- 1
t <- seq(0, 1, by = 0.05)
lambda <- 1 # fsvm & pi.path

# storage
result <- matrix(0, n.sim, 1)

for (iter in 1:n.sim) {
  # tic <- Sys.time()
  seed <- iter
  
  # Data generation
  set.seed(iter)
  data <- gp.1dim(n, beta, t, seed)
  x <- data$x
  y <- data$y
  
  print(iter)
  obj <- fsvm.1dim(y, x, t, L = 10, lambda, rho = 1, weight = rep(1, n))
  K <- obj$K
  
  # pi path
  
  obj_pi <- fsvm.pi.path(lambda, x, y, K)
  
  pi <- obj_pi$pi
  alpha <- obj_pi$alpha
  alpha0 <- matrix(obj_pi$alpha0, dim(alpha)[1], dim(alpha)[2], byrow=T)
  
  new.gx <- K %*% (alpha * y)
  new.fx <- (alpha0 + new.gx)/lambda

  pi.star <- rep(0, n)

  # pi.star
  for (i in 1:n) {
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
  boxplot(pi.star[y == 1], pi.star[y != 1], xlab=iter)
  # dev.off()
  
  # Numerical criteria
  Deviance <- sum(y*log(pi.star))
  Entropy <- -sum(pi.star*log(pi.star)) # 둘중 뭐가 맞는거지?????????????????????????????????
  
  # Deviance
  
  # results  
  result[iter,] <- c(Deviance)
  
  # toc <- Sys.time()
  # print(toc - tic)
}

# save results
write(result, "eu/result/result.txt")


####================================= <Simulation 2> =====================================
# set up
n.sim <- 100
n <- 100
beta <- 1
t <- seq(0, 1, by = 0.05)
lambda <- 1 # fsvm & pi.path
sd <- 0.3

# storage for result
result1 <- matrix(0, n.sim, 1)
Eresult <- matrix(0, n.sim, 1)


for (iter in 1:n.sim) {
  # tic <- Sys.time()
  seed <- iter
  
  # Data generation
  set.seed(iter)
  # create sine function
  n <- 50
  value <- matrix(sin(2*3.14*t), n, length(t), byrow=T)
  
  value[1:25,] <- value[1:25,] + matrix(rnorm(25*length(t), mean=0, sd=sd), 25, length(t), byrow=T)
  value[26:50,] <- value[26:50,] + matrix(rnorm(25*length(t), mean=1, sd=sd), 25, length(t), byrow=T)
  
  value1 <- data.frame(t(value))
  x <- x.list <- as.list(value1)
  y <- c(rep(-1,25),rep(1,25))
  
  print(iter)
  obj <- fsvm.1dim(y, x, t, L = 10, lambda, rho = 1, weight = rep(1, n))
  K <- obj$K
  
  # pi path
  obj_pi <- fsvm.pi.path(lambda, x, y, K)
  
  pi <- obj_pi$pi
  alpha <- obj_pi$alpha
  alpha0 <- matrix(obj_pi$alpha0, dim(alpha)[1], dim(alpha)[2], byrow=T)
  
  new.gx <- K %*% (alpha * y)
  new.fx <- (alpha0 + new.gx)/lambda

  pi.star <- rep(0, n)
  
  # pi.star
  for (i in 1:n) {
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
  boxplot(pi.star[y == 1], pi.star[y != 1], xlab=iter)
  # dev.off()
  
  # Numerical criteria
  Deviance <- sum(y*log(pi.star))
  Entropy <- -sum(pi.star*log(pi.star)) # 둘중 뭐가 맞는거지?????????????????????????????????
  
  # Deviance
  
  # results  
  result1[iter,] <- c(Deviance)
  Eresult[iter,] <- c(Entropy)
  
  # toc <- Sys.time()
  # print(toc - tic)
}

result1
Eresult

# save results
write(result1, "eu/result/result1.txt")







