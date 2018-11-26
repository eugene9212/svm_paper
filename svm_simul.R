#### Final Simulation Code ####
#### for Aug. 24th meeting ####

rm(list = ls())

#### load packages & R code ####
library(mvtnorm)
library(fda)

setwd('C:/Users/eugene/Desktop/SVM/shared/R code/')
source('eu/data_gen/gp.1dim.I.R')     # GP with identity covariance (beta)
source('eu/data_gen/gp.1dim.cov1.R')  # GP with trivial covariance  (beta)
source('eu/data_gen/gp.1dim.sc.R')  # GP with trivial covariance  (beta)
source('eu/data_gen/gp.1dim.AR.R')    # GP with trivial covariance  (beta)
source('eu/data_gen/sc.R')            # non-linear with cross       (no beta)   
source('eu/data_gen/ss.R')            # non-linear with parallel    (beta)
source('eu/data_gen/linear.cross.R')  # linear with cross           (no beta)
source('eu/data_gen/linear.par.R')    # linear with parallel        (beta)

source('eu/fsvm.prob.R')
source('eu/predict.fsvm.prob.R')

source('fn/fsvm.pi.path.R')
source('fn/fsvm.sub.pi.path.R')

dyn.load("C:/Users/eugene/Desktop/SVM_R/shared/R code/KernSurf/temp/wsvmqp.dll")
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('C:/Users/eugene/Desktop/SVM_R/shared/R code/KernSurf/R')


# set up
n.sim <- 50
t <- seq(0, 1, by = 0.05)
# lambda <- 0.5 # fsvm & pi.path
L <- 10
beta <- 0.5

# storage
CRE.result<-matrix(0, n.sim, 1)
pi.result <- as.list(1:n.sim)
lambda.result<-matrix(0, n.sim, 1)
## for test.y
ans <- as.list(1:n.sim)

####========================= Simluation ==================================####
for (iter in 1:n.sim) {
  # iter<-47
  n.train <- 50
  n.test <- 30
  n <- n.train + n.test
  
  # Data generation (6 methods)
  set.seed(iter)
  data <- gp.1dim.I(n, beta, t = t, seed = iter)
  # data <- gp.1dim.cov1(n, beta, t = t, seed = iter)
  # data <- gp.1dim.sc(n, t = t, seed = iter)
  # data <- gp.1dim.AR(n, beta, p = 5, rho = 0.5, t = t, seed = iter)
  # data <- sc(n, t = t, seed = iter)
  # data <- ss(n, beta, t = t, seed = iter)
  # data <- linear.cross(n, t = t, seed = iter)
  # data <- linear.par(n, beta, t = t, seed = iter)
  
  id <- sample(1:n, n.train)
  
  train.x <- data$x[id]
  train.y <- data$y[id]
  test.x <- data$x[-id]
  test.y <- data$y[-id]
  
  print(iter)
  
  ####=======================     train data      ===============================####
  # calculate the pi path
  obj <- fsvm.prob(train.x, train.y, t, L)
  lambda.result[iter,] <- obj$lambda
  
  ####======================= predict probability ===============================####
  obj2 <- predict.fsvm.prob(obj, test.x)
  prob <- obj2$prob

  # Boxplot of pi.star 
  boxplot(prob[test.y == 1], prob[test.y == -1], xlab=paste("SVM",iter), ylim=c(0,1))
  
  # Criteria
  # delta = 1.0e-8
  CRE <- -1/length(test.y)*(sum(1/2*(1+test.y)*log(prob) + 1/2*(1-test.y)*log(1-prob)))
  
  # pi.star
  pi.result[[iter]] <- prob
  
  # Answer
  ans[[iter]] <- test.y
  
  # results  
  CRE.result[iter,] <- c(CRE)
}

# Cross Entropy 
mean(CRE.result)
# median(CRE.result)
# sum(CRE.result)

# Accuracy
pi.result1 <- as.list(1:n.sim)
for (k in 1:n.sim){
  pi.result1[[k]] <- ifelse(pi.result[[k]]<0.5,-1,1)
}

svm <- matrix(0, ncol = n.sim)
for(k in 1:n.sim){
  svm[,k] <- sum(pi.result1[[k]] == ans[[k]])
}

# total number each sim, 30 test.n
mean(svm/30)

paste0("The accuracy of FSVM is ", round(svm/(n.sim*n.test),digits = 3))





