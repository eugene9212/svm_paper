#### Functional Logistic code (binary) ####

rm(list = ls())

#### load packages & R code ####
library(mvtnorm)
library(fda)
library(fda.usc)

setwd('C:/Users/eugene/Desktop/SVM/shared/R code/')
source('eu/data_gen/binary/gp.1dim.I.R')     # GP with identity covariance (beta)
source('eu/data_gen/binary/gp.1dim.cov1.R')  # GP with trivial covariance  (beta)
source('eu/data_gen/binary/gp.1dim.sc.R')  # GP with trivial covariance  (beta)
source('eu/data_gen/binary/gp.1dim.AR.R')    # GP with trivial covariance  (beta)
source('eu/data_gen/binary/gp.1dim.ss.R')  # GP with trivial covariance  (beta)
source('eu/data_gen/binary/linear.cross.R')  # linear with cross           (no beta)
source('eu/data_gen/binary/linear.par.R')    # linear with parallel        (beta)

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
# set up
n.sim <- 50
t <- seq(0, 1, by = 0.03)
rangeval <- quantile(t, c(0,1))
L <- 10
beta <- 1
error <- 0.5

# storage
CRE.fl.result<-matrix(0, n.sim, 1)
pi.fl.result <- as.list(1:n.sim)

####========================= Simluation ==================================####
for (iter in 1:n.sim) {
  # iter<-1
  n.train <- 50
  n.test <- 40
  n <- n.train + n.test
  
  # Data generation (6 methods)
  set.seed(iter)
  data <- gp.1dim.I(n, error, beta, t = t, seed = iter)
  # data <- gp.1dim.cov1(n, error, beta, t = t, seed = iter)
  # data <- gp.1dim.sc(n, error, t = t, seed = iter)
  # data <- gp.1dim.ss(n, error, beta, t = t, seed = iter)
  # data <- gp.1dim.AR(n, error, beta, p = 5, rho = 0.5, t = t, seed = iter)
  # data <- linear.cross(n, error, t = t, seed = iter)
  # data <- linear.par(n, error, beta, t = t, seed = iter)
  
  id <- sample(1:n, n.train)
  
  train.x <- data$x[id]
  train.y <- data$y[id]
  test.x <- data$x[-id]
  test.y <- data$y[-id]
  
  print(iter)
  
  #### Transform train.x and train.y
  ## train.y
  train.fy <- ifelse(train.y == 1, 1, 0)
  dataf <- as.data.frame(train.fy) # transform y as data.frame structure
  
  # train.x -> train.fx(fdata)
  train.f.x.matrix <- matrix(unlist(train.x), nrow = length(train.x), byrow = T)
  x <- train.fx <- fdata(train.f.x.matrix,argvals=t,rangeval=range(t))
  
  #### FDA
  nbasis.x=L # create basis used for fdata or fd covariates.
  nbasis.b=L  # create basis used for beta parameter estimation.
  
  # create basis
  basis1=create.bspline.basis(rangeval=range(t),nbasis=nbasis.x)
  basis2=create.bspline.basis(rangeval=range(t),nbasis=nbasis.b)
  
  # Create basis n ldata before fitting the model
  basis.x=list("x"=basis1) # has to be the same name
  basis.b=list("x"=basis2)
  
  # as.factor
  # dataf <- dataf$train.fy
  ldata=list("df"=dataf,"x"=x)
  
  # formula
  f=train.fy~x 
  
  # Fit the model
  fl.fit <- fregre.glm(f,familiy=binomial(), data=ldata, basis.x=basis.x, basis.b=basis.b, control =list(maxit=1000))
  summary(fl.fit)
  fl.fit$fitted.values
  ####=======================     test data      ===============================####
  # test.x -> test.fx(fdata)
  test.fx.matrix <- matrix(unlist(test.x), nrow = length(test.x), byrow = T)
  test.fx <- fdata(test.fx.matrix, argvals=t, rangeval=range(t))
  
  # create newldata
  newldata <- list("x"=test.fx)
  
  # predict
  pred.glm <- predict.fregre.glm(fl.fit, newx=newldata)
  prob <- pred.glm
  prob <- ifelse(prob < 0,0,prob)
  prob <- ifelse(prob > 1,1,prob)
  
  # Box Plot
  boxplot(prob[test.y == 1], prob[test.y == -1], xlab=paste("logit",iter), ylim=c(0,1))
  
  # Criteria
  CRE <- -1/length(test.Y)*(sum(test.Y*log(prob) + (1-test.Y)*log(1-prob)))
  
  pi.result[[iter]] <- prob
  
  # results  
  CRE.result[iter,] <- c(CRE) 
}

mean(CRE.result)
# median(CRE.result)
# sum(CRE.result)



