#### Final Simulation Code ####
#### for Aug. 24th meeting ####

rm(list = ls())

#### load packages & R code ####
library(mvtnorm)
library(fda)
library(fda.usc)

setwd('C:/Users/eugene/Desktop/SVM/shared/R code/')
source('eu/data_gen/binary/gp.1dim.sc.R')  # GP with trivial covariance  (beta)
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
n.sim <- 100
t <- seq(0, 1, by = 0.05)
L <- 10
error <- 1
beta <- 1
# ridge <- 1e-3

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
  # iter<-92
  n.train <- 500
  n.test <- 1000
  n <- n.train + n.test
  
  # Data generation (4 methods with cov=Identity)
  set.seed(iter)
  # with covariance = Identity
  # data <- gp.1dim.ss(n, error, beta, cov = "I", t = t, seed = iter)
  # data <- gp.1dim.sc(n, error, cov = "I", t = t, seed = iter)
  # data <- linear.cross(n, error, cov = "I", t = t, seed = iter)
  # data <- linear.par(n, error, beta, cov = "I", t = t, seed = iter)
    # with covariance = AR
  data <- gp.1dim.ss(n, error, beta, cov = "AR", t = t, seed = iter)
  # data <- gp.1dim.sc(n, error, cov = "AR", t = t, seed = iter)
  # data <- linear.cross(n, error, cov = "AR", t = t, seed = iter)
  # data <- linear.par(n, error, beta, cov = "AR", t = t, seed = iter)
  # with covariance = CS
  # data <- gp.1dim.ss(n, error, beta, cov = "CS", t = t, seed = iter)
  # data <- gp.1dim.sc(n, error, cov = "CS", t = t, seed = iter)
  # data <- linear.cross(n, error, cov = "CS", t = t, seed = iter)
  # data <- linear.par(n, error, beta, cov = "CS", t = t, seed = iter)
  
  id <- sample(1:n, n.train)
  
  train.x <- data$x[id]
  train.y <- data$y[id]
  test.x <- data$x[-id]
  test.y <- data$y[-id]
  
  ans.p[[iter]] <- data$true.p[-id]
  
  print(iter)
  
  #========================================== Functional SVM ====================####
  ####=======================     train data 
  svm.obj <- tryCatch({
    fsvm.prob(train.x, train.y, t, L)
  }, error = function(e) {
    fsvm.prob(train.x, train.y, t, L, ridge=1e-3)
  })
  # Error in solve.default(Kstar, al)
  # svm.obj <- fsvm.prob(train.x, train.y, t, L, ridge) # calculate the pi path
  
  ####=======================     test data
  svm.obj2 <- predict.fsvm.prob(svm.obj, test.x)
  svm.prob <- svm.obj2$prob
  # cbind(test.y,svm.prob)
  
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
  prob <- ifelse(prob < 0,1e-10,prob)
  prob <- ifelse(prob > 1,1-1e-10,prob)

  # Predicted Probability of Funtional logistic
  pi.fl.result[[iter]] <- c(prob)
  # cbind(test.y, prob)
  
  # Criteria (CRE)
  delta <- 1e-10
  CRE.fl <- -1/length(test.y)*(sum(1/2*(1+test.y)*log(prob+delta) + 1/2*(1-test.y)*log(1-prob+delta)))
  CRE.fl.result[iter,] <- CRE.fl
  
  # Predicted Probability of FSVM
  pi.svm.result[[iter]] <- c(svm.prob)
  
  # Criteria (CRE)
  CRE.svm <- -1/length(test.y)*(sum(1/2*(1+test.y)*log(svm.prob+delta) + 1/2*(1-test.y)*log(1-svm.prob+delta)))
  CRE.svm.result[iter,] <- CRE.svm
  
  # Box Plot
  # boxplot(svm.prob[test.y == 1], svm.prob[test.y == -1], xlab=paste("svm",iter), ylim=c(0,1))
  # boxplot(prob[test.y == 1], prob[test.y == -1], xlab=paste("logit",iter), ylim=c(0,1))

  # store the answers
  ans[[iter]] <- test.y
}

# Check warnings
# summary(warnings())

# Critereon (1) Cross Entropy ======================================================
fl.cre <- round(mean(CRE.fl.result),digits = 3)
svm.cre <- round(mean(CRE.svm.result),digits = 3)

# Critereon (2) Accuracy ===========================================================
svm <- matrix(0, ncol = n.sim)
flog <- matrix(0, ncol = n.sim)
for(k in 1:n.sim){
  flog[,k] <- sum(ifelse(pi.fl.result[[k]]<=0.5, -1, 1) == ans[[k]])
  svm[,k] <- sum(ifelse(pi.svm.result[[k]]<=0.5, -1, 1) == ans[[k]])
}
svm.acc <- round(mean(svm/n.test), digits = 3)
fl.acc <- round(mean(flog/n.test), digits = 3)

# Critereon (3) Distance btw true p & hat p =========================================
predict.p.svm <- pi.svm.result
predict.p.fl <- pi.fl.result

for(i in 1:n.sim){
  # i<-1
  class.fl <- c(); idx.fl <- c(); class.svm <- c(); idx.svm <- c()
  
  class.fl <- ifelse(pi.fl.result[[i]]<=0.5, -1, 1)
  idx.fl <- which(class.fl==-1)
  
  class.svm <- ifelse(pi.svm.result[[i]]<=0.5, -1, 1)
  idx.svm <- which(class.svm==-1)
  
  predict.p.fl[[i]][idx.fl] <- 1-pi.fl.result[[i]][idx.fl]
  predict.p.svm[[i]][idx.svm] <- 1-pi.svm.result[[i]][idx.svm]
}

# calculate the difference
diff.svm <- rep(0,n.sim)
diff.fl <- rep(0,n.sim)
w.diff.svm <- rep(0,n.sim)
w.diff.fl <- rep(0,n.sim)

ans.p1 <- ans.p
# ans.p[[5]]
# i<-5
for (i in 1:n.sim){
  idx1 <- which(c(ans.p1[[i]]-1)==0)
  ans.p1[[i]][idx1] <- c(1-delta)
  
  idx2 <- which(c(ans.p1[[i]]-0)==0)
  ans.p1[[i]][idx2] <- c(delta)
  
  weight <- sqrt(ans.p1[[i]]*(1-ans.p1[[i]]))
  # print(i)
  # print(weight)
  w.diff.svm[i] <- mean(weight*abs(ans.p[[i]] - predict.p.svm[[i]]))
  w.diff.fl[i] <- mean(weight*abs(ans.p[[i]] - predict.p.fl[[i]]))
  
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
paste0("--------------Criterieon (3.1) p diff --------------")
paste0("mean(Diffence) btw true p and svm.predicted p is ", round(mean(diff.svm), digits = 3))
paste0("mean(Diffence) btw true p and fl.predicted p is ", round(mean(diff.fl), digits = 3))
paste0("--------------Criterieon (3.2) p diff (weighted) --------------")
paste0("mean(Diffence) btw true p and svm.predicted p is ", round(mean(w.diff.svm), digits = 3))
paste0("mean(Diffence) btw true p and fl.predicted p is ", round(mean(w.diff.fl), digits = 3))
alarm()