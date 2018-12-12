#### Final Simulation Code ####
#### for Aug. 24th meeting ####

rm(list = ls())

#### load packages & R code ####
library(mvtnorm)
library(fda)

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
sourceDir('C:/Users/eugene/Desktop/SVM_R/shared/R code/KernSurf/R')


# set up
n.sim <- 30
t <- seq(0, 1, by = 0.05)
L <- 10
error <- 0.5
beta <- 1

# storage
## for test.y
ans <- as.list(1:n.sim)
ans.p <- as.list(1:n.sim)
## for Functional logistic
CRE.fl.result<-matrix(0, n.sim, 1)
pi.fl.result <- as.list(1:n.sim)
## for SVM
CRE.svm.result<-matrix(0, n.sim, 1)
pi.svm.result <- as.list(1:n.sim)

####========================= Simluation ==================================####
for (iter in 1:n.sim) {
  # iter<-47
  n.train <- 60
  n.test <- 40
  n <- n.train + n.test
  
  # Data generation (6 methods)
  set.seed(iter)
  data <- gp.1dim.I(n, error, beta, t = t, seed = iter)
  # data <- gp.1dim.cov1(n, error, beta, t = t, seed = iter)
  # data <- gp.1dim.sc(n, error, t = t, seed = iter)
  # data <- gp.1dim.ss(n, error, t = t, seed = iter)
  # data <- gp.1dim.AR(n, error, beta, p = 5, rho = 0.5, t = t, seed = iter)
  # data <- linear.cross(n, error, t = t, seed = iter)
  # data <- linear.par(n, error, beta, t = t, seed = iter)
  
  id <- sample(1:n, n.train)
  
  train.x <- data$x[id]
  train.y <- data$y[id]
  test.x <- data$x[-id]
  test.y <- data$y[-id]
  
  print(iter)
  
  #========================================== Functional SVM ====================####
  ####=======================     train data 
  # calculate the pi path
  svm.obj <- fsvm.prob(train.x, train.y, t, L)
  
  ####=======================     test data
  svm.obj2 <- predict.fsvm.prob(svm.obj, test.x)
  svm.prob <- svm.obj2$prob
  
  #========================================== Functional logistic ===============####
  #### Transform train.x and train.y
  train.fy <- ifelse(train.y == 1, 1, 0)
  train.fy <- data.frame(train.fy)
  train.f.x.matrix <- matrix(unlist(train.x), nrow = length(train.x), byrow = T)
  train.fx <- fdata(train.f.x.matrix,argvals=t,rangeval=range(t))
  
  #### FDA
  nbasis.x=L; nbasis.b=L
  
  # create basis
  basis1=create.bspline.basis(rangeval=range(t),nbasis=nbasis.x)
  basis2=create.bspline.basis(rangeval=range(t),nbasis=nbasis.b)
  # Create basis n ldata before fitting the model
  basis.x=list("x"=basis1) ; basis.b=list("x"=basis2)
  
  # formula
  f=train.fy~x 
  # as.factor
  train.fy <- train.fy$train.fy
  ldata=list("df"=train.fy,"x"=train.fx)
  
  ####=======================     train data 
  # Fit the model
  fl.fit <- fregre.glm(f,familiy=binomial(link = "logit"), data=ldata, basis.x=basis.x, basis.b=basis.b, control =list(maxit=1000))
  
  ####=======================     test data  
  # test.x -> test.fx(fdata)
  test.fx.matrix <- matrix(unlist(test.x), nrow = length(test.x), byrow = T)
  test.fx <- fdata(test.fx.matrix, argvals=t, rangeval=range(t))
  
  # create newldata
  newldata <- list("x"=test.fx)
  
  # predict
  prob <- predict.fregre.glm(fl.fit, newldata)
  fl.prob <- exp(prob)/(1+exp(prob))
  
  # Predicted Probability of Funtional logistic
  pi.fl.result[[iter]] <- fl.prob
  
  # Criteria (CRE)
  CRE.fl <- -1/length(test.y)*(sum(1/2*(1+test.y)*log(fl.prob) + 1/2*(1-test.y)*log(1-fl.prob)))
  CRE.fl.result[iter,] <- CRE.fl
  
  # Predicted Probability of FSVM
  pi.svm.result[[iter]] <- svm.prob
  
  # Criteria (CRE)
  CRE.svm <- -1/length(test.y)*(sum(1/2*(1+test.y)*log(svm.prob) + 1/2*(1-test.y)*log(1-svm.prob)))
  CRE.svm.result[iter,] <- CRE.svm
  
  # store the answers
  ans[[iter]] <- test.y
  ans.p[[iter]] <- data$true.p[-id]
  
}

# Check warnings
summary(warnings())

# Critereon (1) Cross Entropy
fl.cre<-round(mean(CRE.fl.result),digits = 3)
svm.cre<-round(mean(CRE.svm.result),digits = 3)

# Critereon (2) Accuracy
svm <- matrix(0, ncol = n.sim)
flog <- matrix(0, ncol = n.sim)
for(k in 1:n.sim){
  svm[,k] <- sum(apply(pi.fl.result[[k]], 1, which.max) == ans[[k]])
  flog[,k] <- sum(apply(pi.svm.result[[k]], 1, which.max) == ans[[k]])
}

# total number : 50 * 30 = 1500
svm.acc <- round(mean(svm/n.test), digits = 3)
fl.acc <- round(mean(flog/n.test), digits = 3)

# Critereon (3) Distance btw true p & hat p
predict.p.svm <- as.list(1:n.sim)
predict.p.fl <- as.list(1:n.sim)

for(i in 1:n.sim){
  idx <- ans[[i]]
  for(j in 1:n.test){
    a <- pi.svm.result[[i]][j,]
    b <- pi.fl.result[[i]][j,]
    predict.p.svm[[i]][j] <- a[idx[j]] 
    predict.p.fl[[i]][j] <- b[idx[j]] 
  }
}

# calculate the difference
diff.svm <- rep(0,n.sim)
diff.fl <- rep(0,n.sim)

for (i in 1:n.sim){
  diff.svm[i] <- mean(abs(ans.p[[i]] - predict.p.svm[[i]]))
  diff.fl[i] <- mean(abs(ans.p[[i]] - predict.p.fl[[i]]))
}

# print
paste0("--------------Simulation Result with Error ",error,"--------------")
paste0("--------------Criterieon (1) CRE --------------")
paste0("CRE of FSVM is ", round(svm.cre, digits = 3))
paste0("CRE of Flogistic is ", round(fl.cre, digits = 3))
paste0("--------------Criterieon (2) Accuracy --------------")
paste0("Accuracy of FSVM is ", round(svm.acc, digits = 3))
paste0("Accuracy of Flogistic is ", round(fl.acc, digits = 3))
paste0("--------------Criterieon (3) p diff --------------")
paste0("mean(Diffence) btw true p and svm.predicted p is ", round(mean(diff.svm), digits = 3))
paste0("mean(Diffence) btw true p and fl.predicted p is ", round(mean(diff.fl), digits = 3))




