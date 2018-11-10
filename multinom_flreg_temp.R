#### Class Simulation Code   ####
##Multinomial Functional logistic   ####

rm(list = ls())

#### load packages & R code ####
library(mvtnorm)
library(fda)
library(fda.usc)

setwd('C:/Users/eugene/Desktop/SVM/shared/R code/')
source('eu/data_gen/multiclass/class3/gp.I.crss.linear.3.error.R')
source('eu/data_gen/multiclass/class3/gp.I.linear.K.error.R')
source('eu/data_gen/multiclass/class3/gp.I.nonlinear.3.error.R')

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
rangeval <- quantile(t, c(0,1))
L <- 10
beta <- 1
error <- 0.1
K <- 3

# storage
CRE.result<-matrix(0, n.sim, 1)
pi.result <- as.list(1:n.sim)
# iter.result<-matrix(0, n.sim, 1)

####========================= Simluation ==================================####
for (iter in 1:n.sim) {
  iter<-47
  n.train <- 60
  n.test <- 30
  n <- n.train + n.test
  
  # Data generation (3 methods)
  set.seed(iter)
  data <- gp.I.crss.linear.3.error(n, error, t, seed = iter)
  # data <- gp.I.linear.K.error(n, error, beta, K, t, seed = iter)
  # data <- gp.I.nonlinear.3.error(n, error, t, seed = iter)
  
  id <- sample(1:n, n.train)
  
  train.x <- data$x[id]
  train.y <- data$y[id]
  test.x <- data$x[-id]
  test.y <- data$y[-id]
  
  print(iter)
  
  # as.factor(y)
  train.y <- as.factor(train.y)
  test.y <- as.factor(test.y)
  
  #### Transform train.x and train.y
  train.f.x <- train.x
  train.f.x.matrix <- matrix(unlist(train.f.x), nrow = length(train.f.x), byrow = T)
  D <- list("data" = train.f.x.matrix, "argvals" = t, rangeval = rangeval)
  attr(D, "class") <- "fdata"
  
  train.f.y <- train.y
  train.f.y <- as.data.frame(train.f.y)
  
  basis.obj <- create.bspline.basis(rangeval, 10) # L : number of basis
  basis.x <- list("x"=basis.obj)
  f <- train.f.y ~ x
  ldata <- list("df"=train.f.y,"x"=D)
  
  ####=======================     train data      ===============================####
  res <- fregre.gsam(f, family=multinomial(), data=ldata, basis.x=basis.x,maxit = 1000)
  summary(res)
  ?multinom
  
  # Predict
  test.X <- test.x
  test.x.matrix <- matrix(unlist(test.X), nrow = length(test.X), byrow = T)
  test.D <- list("data" = test.x.matrix, "argvals" = t, rangeval = rangeval)
  attr(test.D, "class") <- "fdata"
  
  test.Y <- test.y; index <- test.Y<0; test.Y[index] = 0
  test.dataf <- as.data.frame(test.Y)
  
  f <- test.Y ~ x
  newldata <- list("df"=test.dataf,"x"=test.D)
  
  # predict(res, newldata, type = "response")
  pred.glm <- predict.fregre.glm(res, newldata)
  prob <- pred.glm
  
  boxplot(prob[test.Y == 1], prob[test.Y == 0], xlab=paste("logit",iter), ylim=c(0,1))
  
  # Criteria
  CRE <- -1/length(test.Y)*(sum(test.Y*log(prob) + (1-test.Y)*log(1-prob)))
  
  pi.result[[iter]] <- prob
  
  # results  
  CRE.result[iter,] <- c(CRE) 
}

mean(CRE.result)
# median(CRE.result)
# sum(CRE.result)



