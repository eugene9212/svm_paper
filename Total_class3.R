#### Class Simulation Code   ####
####   Functional logistic   ####

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
n.sim <- 100
t <- seq(0, 1, by = 0.05)
rangeval <- quantile(t, c(0,1))
L <- 10
beta <- 1
error <- 1
K <- 3
max.iter <- 100

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
  # iter<-1
  n.train <- 60
  n.test <- 30
  n <- n.train + n.test
  
  # Data generation (3 methods)
  set.seed(iter)
  # data <- gp.I.linear.K.error(n, error, beta, K, t, seed = iter)
  data <- gp.I.crss.linear.3.error(n, error, t, seed = iter)
  # data <- gp.I.nonlinear.3.error(n, error, t, seed = iter)
  
  id <- sample(1:n, n.train)
  
  train.x <- data$x[id]
  train.y <- data$y[id]
  test.x <- data$x[-id]
  test.y <- data$y[-id]
  
  true.p <- data$true.p[-id]
  
  print(iter)
  
  # Divide train set pairwise
  trn.idx12 <- which(train.y==1 | train.y==2)
  trn.idx13 <- which(train.y==1 | train.y==3)
  trn.idx23 <- which(train.y==2 | train.y==3)
  
  # train.x
  train.x.12 <- train.x[trn.idx12]
  train.x.13 <- train.x[trn.idx13]
  train.x.23 <- train.x[trn.idx23]
  
  # train.y
  train.y.12 <- train.y[trn.idx12]
  train.y.13 <- train.y[trn.idx13]
  train.y.23 <- train.y[trn.idx23]
  
  # Change y value into -1 and 1
  train.y.12 <- ifelse(train.y.12 == 2, 1, -1)
  train.y.13 <- ifelse(train.y.13 == 3, 1, -1)
  train.y.23 <- ifelse(train.y.23 == 3, 1, -1)
  
  #### Transform train.x and train.y
  # train.x.12 -> train.fx12(fdata)
  train.f.x.12.matrix <- matrix(unlist(train.x.12), nrow = length(train.x.12), byrow = T)
  train.fx12 <- fdata(train.f.x.12.matrix,argvals=t,rangeval=range(t))
  # train.x.13 -> train.fx13(fdata)
  train.f.x.13.matrix <- matrix(unlist(train.x.13), nrow = length(train.x.13), byrow = T)
  train.fx13 <- fdata(train.f.x.13.matrix,argvals=t,rangeval=range(t))
  # train.x.23 -> train.fx23(fdata)
  train.f.x.23.matrix <- matrix(unlist(train.x.23), nrow = length(train.x.23), byrow = T)
  train.fx23 <- fdata(train.f.x.23.matrix,argvals=t,rangeval=range(t))
  
  # train.y.12 -> train.fy12 (0,1)
  train.fy12 <- train.y.12
  index12 <- train.y.12 < 0
  train.fy12[index12] = 0
  train.fy12 <- as.data.frame(train.fy12) # transform y as data.frame structure
  # train.y.13 -> train.fy13 (0,1)
  train.fy13 <- train.y.13
  index13 <- train.y.13 < 0
  train.fy13[index13] = 0
  train.fy13 <- as.data.frame(train.fy13) # transform y as data.frame structure
  # train.y.23 -> train.fy23 (0,1)
  train.fy23 <- train.y.23
  index23 <- train.y.23 < 0
  train.fy23[index23] = 0
  train.fy23 <- as.data.frame(train.fy23) # transform y as data.frame structure
  
  #### FDA
  nbasis.x=L # create basis used for fdata or fd covariates.
  nbasis.b=L  # create basis used for beta parameter estimation.
  
  # create basis
  basis1=create.bspline.basis(rangeval=range(t),nbasis=nbasis.x)
  basis2=create.bspline.basis(rangeval=range(t),nbasis=nbasis.b)
  
  # formula
  f12=train.fy12~x 
  f13=train.fy13~x 
  f23=train.fy23~x 
  
  # Create basis n ldata before fitting the model
  basis.x=list("x"=basis1) # has to be the same name
  basis.b=list("x"=basis2)
  
  # as.factor
  train.fy12 <- train.fy12$train.fy12
  train.fy13 <- train.fy13$train.fy13
  train.fy23 <- train.fy23$train.fy23
  ldata12=list("df"=train.fy12,"x"=train.fx12)
  ldata13=list("df"=train.fy13,"x"=train.fx13)
  ldata23=list("df"=train.fy23,"x"=train.fx23)
  
  # Fit the model
  res12=fregre.glm(f12,familiy=binomial(link = "logit"), data=ldata12, basis.x=basis.x, basis.b=basis.b, control =list(maxit=1000))
  res13=fregre.glm(f13,familiy=binomial(link = "logit"), data=ldata13, basis.x=basis.x, basis.b=basis.b, control =list(maxit=1000))
  res23=fregre.glm(f23,familiy=binomial(link = "logit"), data=ldata23, basis.x=basis.x, basis.b=basis.b, control =list(maxit=1000))
  
  ####=======================     test data      ===============================####
  # test.x -> test.fx(fdata)
  test.fx.matrix <- matrix(unlist(test.x), nrow = length(test.x), byrow = T)
  test.fx <- fdata(test.fx.matrix, argvals=t, rangeval=range(t))
  
  # create newldata
  newldata <- list("x"=test.fx)
  
  # predict
  pred.glm12 <- predict.fregre.glm(res12, newldata)
  pred.glm12 <- ifelse(pred.glm12 < 0,1e-10,pred.glm12)
  pred.glm12 <- ifelse(pred.glm12 > 1,1-1e-10,pred.glm12)
  # test.y
  # round(pred.glm12)
  pred.glm13 <- predict.fregre.glm(res13, newldata)
  pred.glm13 <- ifelse(pred.glm13 < 0,1e-10,pred.glm13)
  pred.glm13 <- ifelse(pred.glm13 > 1,1-1e-10,pred.glm13)
  # test.y
  # round(pred.glm13)
  pred.glm23 <- predict.fregre.glm(res23, newldata)
  pred.glm23 <- ifelse(pred.glm23 < 0,1e-10,pred.glm23)
  pred.glm23 <- ifelse(pred.glm23 > 1,1-1e-10,pred.glm23)
  # test.y
  # round(pred.glm23)
  
  # 시뮬레이션 1번=test sample은 n.test개. 따라서 p도 시뮬 한번 당 n.test개 나옴 -------------------#
  # 시뮬당 저장공간 생성
  pi.fl.one.simul <- matrix(0, n.test, K)
  test.y.one.simul <- matrix(test.y, n.test, 1)
  
  # Pairwise Coupling
  for(ii in 1:n.test){
    # print(ii)
    r <- matrix(0, K, K)
    r[lower.tri(r)] <- c(pred.glm12[ii], # r21
                         pred.glm13[ii], # r31
                         pred.glm23[ii]) # r32
    r[upper.tri(r)] <- t(r)[upper.tri(r)]
    one <- matrix(1,K,K)
    one[lower.tri(one, diag=TRUE)] <- 0
    r <- one - r
    r <- abs(r)
    # 1e-323==0
    
    ## Algorithm2.
    ### Create Q matrix
    Q <- matrix(0,K,K)
    for(i in 1:K){
      for(j in 1:K){
        Q[i,j] <- -r[j,i]*r[i,j]
      }
    }
    diag(Q) <- colSums(r^2)
    
    ### (1) Initialize P
    p <- matrix(rep(1/K),K)
    p[K] <- 1-sum(p[-K])
    
    ### (2) Repeat (tt = 1, 2, 3, ..., K, 1, ...)
    tt <- 1
    iter.n <- 1
    
    while(iter.n <= max.iter){
      a <- 1/Q[tt,tt]
      b <- t(p) %*% Q %*% p
      p[tt,] <- a * ( -as.vector(Q[tt,-tt]) %*% p[-tt]  + b)
      
      ## normalize
      p <- p/sum(p)
      
      ## Condition (21) check
      tmp <- Q %*% p
      tmp2 <- matrix(outer(tmp, tmp, "-"), K, K)
      tmp3 <- tmp2[upper.tri(tmp2)]
      # idx <- which(max(p) == p)
      # c <- c(p[idx] - p[-idx])
      
      if (abs(tmp3)[1] < 1e-05 && abs(tmp3)[2] < 1e-05 && abs(tmp3)[3] < 1e-05) break
      # if (abs(tmp3) < 1e-05 && sum(p) == 1) break
      # if (length(unique(sign(tmp))) == 1 && abs(tmp3) < 1e-05 && sum(p) == 1) break
      
      ## re-indexing t
      if (tt == K) {
        tt <- 1
      } else {
        tt <- (tt + 1)
      }
      
      iter.n <- iter.n + 1 # counting the iteration
    }
    
    if (iter.n == max.iter) warning("maximum iteration reached!")
    
    pi.fl.one.simul[ii,] <- p
    
  } #--- End of n.test loop ------------------------------------------------------------------------#
  
  # 시뮬레이션 1번=test sample은 n.test개. 따라서 p도 시뮬 한번 당 n.test개 나옴 -------------------#
  
  # ------------------------------ Simulation on FSVM ----------------------------------------------#
  ####=======================     train data      ===============================####
  # calculate the pi path
  obj12 <- fsvm.prob(train.x.12, train.y.12, t, L)
  obj13 <- fsvm.prob(train.x.13, train.y.13, t, L)
  obj23 <- fsvm.prob(train.x.23, train.y.23, t, L)
  
  ####======================= predict probability ===============================####
  obj212 <- predict.fsvm.prob(obj12, test.x)
  obj213 <- predict.fsvm.prob(obj13, test.x)
  obj223 <- predict.fsvm.prob(obj23, test.x)
  
  # 시뮬레이션 1번=test sample은 n.test개. 따라서 p도 시뮬 한번 당 n.test개 나옴 -------------------#
  # 시뮬당 저장공간 생성
  pi.svm.one.simul <- matrix(0, n.test, K)
  
  # Pairwise Coupling
  for(ii in 1:n.test){
    r <- matrix(0, K, K)
    r[lower.tri(r)] <- c(obj212$prob[ii], # r21
                         obj213$prob[ii], # r31
                         obj223$prob[ii]) # r32
    r[upper.tri(r)] <- t(r)[upper.tri(r)]
    one <- matrix(1,K,K)
    one[lower.tri(one, diag=TRUE)] <- 0
    r <- one - r
    r <- abs(r)
    
    ## Algorithm2.
    ### Create Q matrix
    Q <- matrix(0,K,K)
    for(i in 1:K){
      for(j in 1:K){
        Q[i,j] <- -r[j,i]*r[i,j]
      }
    }
    diag(Q) <- colSums(r^2)
    
    ### (1) Initialize P
    p <- matrix(rep(1/K),K)
    p[K] <- 1-sum(p[-K])

    ### (2) Repeat (tt = 1, 2, 3, ..., K, 1, ...)
    tt <- 1
    iter.n <- 1
    
    while(iter.n <= max.iter){
      a <- 1/Q[tt,tt]
      b <- t(p) %*% Q %*% p
      p[tt,] <- a * ( -as.vector(Q[tt,-tt]) %*% p[-tt]  + b)
      
      ## normalize
      p <- p/sum(p)
      
      ## Condition (21) check
      tmp <- Q %*% p
      tmp2 <- matrix(outer(tmp, tmp, "-"), K, K)
      tmp3 <- tmp2[upper.tri(tmp2)]
      # idx <- which(max(p) == p)
      # c <- c(p[idx] - p[-idx])
      
      if (abs(tmp3)[1] < 1e-05 && abs(tmp3)[2] < 1e-05 && abs(tmp3)[3] < 1e-05) break
      # if (abs(tmp3) < 1e-05 && sum(p) == 1) break
      # if (length(unique(sign(tmp))) == 1 && abs(tmp3) < 1e-05 && sum(p) == 1) break
      
      ## re-indexing t
      if (tt == K) {
        tt <- 1
      } else {
        tt <- (tt + 1)
      }
      
      iter.n <- iter.n + 1 # counting the iteration
    }
    
    if (iter.n == max.iter) warning("maximum iteration reached!")
    
    pi.svm.one.simul[ii,] <- p
    
  } #--- End of n.test loop ------------------------------------------------------------------------#
  
  # 시뮬레이션 1번=test sample은 n.test개. 따라서 p도 시뮬 한번 당 n.test개 나옴 -------------------#
  
  
  # Save results (for flogistic) -------------------------------------------------------------------#
  # CRE 
  s.fl <- c()
  for (k in 1:n.test){
    s.fl <- c(s.fl,-log(pi.fl.one.simul[k,test.y[k]]))
  }
  CRE.fl <- sum(s.fl)
  CRE.fl.result[iter] <- CRE.fl
  
  # pi.star
  pi.fl.result[[iter]] <- pi.fl.one.simul 
  
  # Save results (for FSVM) -----------------------------------------------------------------------#
  # CRE 
  s.svm <- c()
  for (k in 1:n.test){
    s.svm <- c(s.svm,-log(pi.svm.one.simul[k,test.y[k]]))
  }
  CRE.svm <- sum(s.svm)
  CRE.svm.result[iter] <- CRE.svm
  
  # Answer
  ans[[iter]] <- test.y
  ans.p[[iter]] <- true.p
  
  # pi.star
  pi.svm.result[[iter]] <- pi.svm.one.simul 
}

# Check warnings
summary(warnings())

# Critereon (1) Cross Entropy
fl.cre<-round(mean(CRE.fl.result),digits = 3)
svm.cre<-round(mean(CRE.svm.result),digits = 3)
# median(CRE.fl.result)
# median(CRE.svm.result)
# sum(CRE.fl.result)
# sum(CRE.svm.result)

# Critereon (2) Accuracy
svm <- matrix(0, ncol = n.sim)
flog <- matrix(0, ncol = n.sim)
for(k in 1:n.sim){
  flog[,k] <- sum(apply(pi.fl.result[[k]], 1, which.max) == ans[[k]])
  svm[,k] <- sum(apply(pi.svm.result[[k]], 1, which.max) == ans[[k]])
}

# total number : 50 * 30 = 1500
svm.acc <- round(mean(svm/n.test), digits = 3)
fl.acc <- round(mean(flog/n.test), digits = 3)

# print
paste0("--------------Simulation Result with Error ",error,"--------------")
paste0("The CRE of Flogistic is ", fl.cre)
paste0("The CRE of FSVM is ", svm.cre)
paste0("The accuracy of Flogistic is ", fl.acc)
paste0("The accuracy of FSVM is ", svm.acc)

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
paste0("CRE of Flogistic is ", round(fl.cre, digits = 3))
paste0("CRE of FSVM is ", round(svm.cre, digits = 3))
paste0("--------------Criterieon (2) Accuracy --------------")
paste0("Accuracy of Flogistic is ", round(fl.acc, digits = 3))
paste0("Accuracy of FSVM is ", round(svm.acc, digits = 3))
paste0("--------------Criterieon (3) p diff --------------")
# paste0("sum(Diffence) btw true p and svm.predicted p is ", round(sum(diff.svm), digits = 3))
# paste0("sum(Diffence) btw true p and fl.predicted p is ", round(sum(diff.fl), digits = 3))
paste0("mean(Diffence) btw true p and svm.predicted p is ", round(mean(diff.svm), digits = 3))
paste0("mean(Diffence) btw true p and fl.predicted p is ", round(mean(diff.fl), digits = 3))
# paste0("--------------Criterieon (4) KL divergence --------------")
# paste0("KL Divergence with svm.p is ", round(mean(kl.svm), digits = 3))
# paste0("KL Divergence with fl.p is ", round(mean(kl.fl), digits = 3))

