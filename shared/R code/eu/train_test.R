
########################################
## Simulation with train and test set ##
########################################

rm(list = ls())

#### load packages
library(mvtnorm)
library(fda)

#### load R codes
setwd('C:/Users/eugene/Desktop/SVM_R/shared/R code/')
source('eu/fsvm.1dim.R')
source('eu/gp.1dim.sc.R')
source('eu/gp.1dim.R')
source('eu/predict.1dim.fsvm.R')
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

####================================= <Simulation> =====================================
# set up
n.sim <- 100
n <- 100
t <- seq(0, 2*pi, by = 0.05)
lambda <- 1 # fsvm & pi.path

# storage
result <- matrix(0, n.sim, 1)

for (iter in 1:n.sim) {
  # tic <- Sys.time()
  iter <- 1
  
  # Data generation
  set.seed(iter)
  data <- gp.1dim.sc(n,t, seed = iter)
  
  # train and test set
  id1 <- sample(n/2, n/4)
  id2 <- sample((n/2+1):n, n/4)
  id <- c(id1, id2)
  
  train.x <- data$x[id]
  train.y <- data$y[id]
  test.x <- data$x[-id]
  test.y <- data$y[-id]
  
  print(iter)
  obj <- fsvm.1dim(train.y, train.x, t, L = 10, lambda, rho = 1, weight = rep(1, n))
  # K <- obj$K
  
  # -------------------- 여기서부터(꼭!!) --------------------------------------- #
  # predict y with test.x ## 이부분 확인필요!! 
  # new.fx를 구해서 음수 양수로 예측후 pi.path하는게 맞는지????????????????????????????????????????????????
  obj.p <- predict.1dim.fsvm(obj, test.x)
  predict.y <- as.vector(ifelse(obj.p$new.fx < 0, -1, 1))
    
  # pi path
  obj_pi <- fsvm.pi.path(lambda, test.x, predict.y, obj.p$new.K) # obj$k인지/ obj.p$new.k인지 
  # y가 필요. 따라서 predict후 pi.path 쓸 수 밖에!!
  # 여기서 K는 어떤거???????????????????????????????????????????(new.k? train후의 k?)
  
  # predict y
  pi <- obj_pi$pi
  alpha <- obj_pi$alpha
  alpha0 <- matrix(obj_pi$alpha0, dim(alpha)[1], dim(alpha)[2], byrow=T)
  
  new.gx <- obj.p$new.K %*% (alpha * predict.y) # 여기서 K는 어떤거??????????????????????????????????(new.k? train후의 k?)
  new.fx <- (alpha0 + new.gx)/lambda
  # -------------------- 여기까지 ----------------------------------------- # 
  pi.star <- rep(0, length(test.y))
  
  # pi.star
  for (i in 1:length(test.y)) {
    i <- 1
    new.fx[i,]
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
  # dev.off()
  
  # Numerical criteria
  CRE <- -1/length(test.y)*sum(1/2*(1+test.y)*log(pi.star) + 1/2*(1-test.y)*log(1-pi.star))
  
  # results  
  result[iter,] <- c(CRE)
  
  # toc <- Sys.time()
  # print(toc - tic)
}

# save results
write(result, "eu/result/result.txt")


