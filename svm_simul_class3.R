#### Class Simulation Code   ####
#### after Oct. 30th meeting ####

rm(list = ls())

#### load packages & R code ####
library(mvtnorm)
library(fda)

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
L <- 10
beta <- 1
error <- 0.1
K <- 3

# storage
CRE.result<-matrix(0, n.sim, 1)
pi.result <- as.list(1:n.sim)
ans <- as.list(1:n.sim)
# iter.result<-matrix(0, n.sim, 1)

####========================= Simluation ==================================####
for (iter in 1:n.sim) {
  # iter<-3
  # t <- seq(0, 1, by = 0.05)
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
  pi.one.simul <- matrix(0, n.test, K)
  test.y.one.simul <- matrix(test.y, n.test, 1)
  
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
    p <- matrix(rep(1/K),K) # 이렇게 초기화해도돼나Q
    p[K] <- 1-sum(p[-K])
    
    ### (2) Repeat (tt = 1, 2, 3, ..., K, 1, ...)
    tt <- 1
    iter.n <- 1
    
    while(TRUE){
      a <- 1/Q[tt,tt]
      b <- t(p) %*% Q %*% p
      p[tt,] <- a * ( -as.vector(Q[tt,-tt]) %*% p[-tt]  + b)
      
      ## normalize
      p <- p/sum(p)
      
      ## Condition (21) check
      tmp <- Q %*% p
      tmp2 <- matrix(outer(tmp, tmp, "-"), K, K)
      tmp3 <- tmp2[upper.tri(tmp2)]
      idx <- which(max(p) == p)
      c <- c(p[idx] - p[-idx])
      
      if (length(unique(sign(tmp))) == 1 && abs(tmp3) < 1e-05 && sum(p) == 1) break
      
      ## re-indexing t
      if (tt == K) {
        tt <- 1
      } else {
        tt <- (tt + 1)
      }
      
      # iter.n <- iter.n + 1 # counting the iteration
    }
    
    pi.one.simul[ii,] <- p
    
  } #--- End of n.test loop ------------------------------------------------------------------------#
  
  # 시뮬레이션 1번=test sample은 n.test개. 따라서 p도 시뮬 한번 당 n.test개 나옴 -------------------#
    
  # Save results
  # CRE 
  s <- c()
  for (k in 1:n.test){
    s <- c(s,-log(pi.one.simul[k,test.y[k]]))
  }
  CRE <- sum(s)
  CRE.result[iter] <- CRE
  
  # Answer
  ans[[iter]] <- test.y
  
  # pi.star
  pi.result[[iter]] <- pi.one.simul
}

# Performance #
CRE.result
mean(CRE.result)

ans

class.1 <- c()
class.2 <- c()
class.3 <- c()
for(i in 1:n.sim){
  for(j in 1:n.test){
    idx1 <- which(ans[[i]] == 1)
    idx2 <- which(ans[[i]] == 2)
    idx3 <- which(ans[[i]] == 3)
    
    class.1 <- c(class.1, pi.result[[i]][idx1,1])
    class.2 <- c(class.2, pi.result[[i]][idx2,2])
    class.3 <- c(class.3, pi.result[[i]][idx3,3])
  }
}

boxplot(class.1)
boxplot(class.2)
boxplot(class.3)




