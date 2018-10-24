rm(list = ls())

#### load packages & R code ####
library(mvtnorm)
library(fda)

setwd('C:/Users/eugene/Desktop/SVM/shared/R code/')
source('eu/data_gen/multiclass/class.3.R')
source('eu/data_gen/multiclass/class.4.R')

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
t <- seq(0, 1, by = 0.03)
L <- 10
error <- 0.9

n.train <- 90
n.test <- 30
n <- n.train + n.test

# Data generation
data <- class.3(n, error, t = t, seed = 1); K <- 3
# data <- class.4(n, error, t = t, seed = 1); K <- 4

id <- sample(1:n, n.train)

train.x <- data$x[id]
train.y <- data$y[id]
test.x <- data$x[-id]
test.y <- data$y[-id]

# Divide train set pairwise
trn.idx12 <- which(train.y==1 | train.y==2)
trn.idx13 <- which(train.y==1 | train.y==3)
trn.idx23 <- which(train.y==2 | train.y==3)
train.x.12 <- train.x[trn.idx12]
train.x.13 <- train.x[trn.idx13]
train.x.23 <- train.x[trn.idx23]
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

# Save Data
save(obj212, obj213, obj223, file="C:/Users/eugene/Desktop/class3.RData")

####============================= Pairwise Calculation =====================================####
###=============================== Algorithm 1. ============================================####
# create pairwise n matrix
n <- matrix(0, K, K) # n
n[upper.tri(n)] <- c(length(train.y.12), # n12
                     length(train.y.13), # n13
                     length(train.y.23)) # n23
n[lower.tri(n)] <- t(n)[lower.tri(n)]
diag(n) <- rep(0,K)

result <- as.list(1:n.test)
r <- matrix(0, K, K)

for (ii in 1:n.test){
  # ii <- 1
  r <- matrix(0, K, K)
  r[lower.tri(r)] <- c(obj212$prob[ii], # r21
                       obj213$prob[ii], # r31
                       obj223$prob[ii]) # r32
  r[upper.tri(r)] <- 1-c(obj212$prob[ii], # r12
                         obj213$prob[ii], # r13
                         obj223$prob[ii]) # r23

  #### Pairwise Calculation ####
  ### Algorithm 1. ####
  ## 1. Initialize pi & corresponding muij
  p <- runif(K,0,1) # pi
  tmp <- outer(p,p,"+")
  mu <- p/tmp # mu
  diag(mu) <- rep(0,K)
  
  check <- c()
  
  ## 2. Repeat (i = 1,2, ... ,k,1, ...)
  i <- 1
  while (TRUE){
    alpha <- sum(n[i,] * r[i,]) / sum(n[i,] * mu[i,])
    mu[i,] <- alpha*mu[i,]/(alpha*mu[i,]+mu[,i])
    mu[,i] <- 1 - mu[i,]
    print(mu)
    diag(mu) <- rep(0,K)
    
    check <- c(check, alpha)

  ## 3. Convergence check
    if (i >= K && abs(1 - check[(i-K+1) : i]) < 1e-08) break
    if (i == K) i <- 1
    else i <- (i + 1)
    
  ## warning
    if (length(check) > 1000) message("Warning : Iteration is over 1000.")
  }
  
  p <- p*alpha
  p <- p/sum(p)
  
  result[[ii]] <-  p
}


####============================= Pairwise Calculation =====================================####
###=============================== Algorithm 2. ============================================####
# load Data
load(file="class3.RData")
K <- 3
n <- 120
n12 <- 60 # length(train.y.12)
n13 <- 60 # length(train.y.12)
n23 <- 60 # length(train.y.12)
n.test <- 30

# create pairwise n matrix
n <- matrix(0, K, K) # n
n[upper.tri(n)] <- c(n12,n13,n23) # n23
n[lower.tri(n)] <- t(n)[lower.tri(n)]
diag(n) <- rep(0,K)

result <- as.list(1:n.test)
r <- matrix(0, K, K)

# Parellel computing
working1<-function(ii){
  r <- matrix(0, K, K)
  r[lower.tri(r)] <- c(obj212$prob[ii], # r21
                       obj213$prob[ii], # r31
                       obj223$prob[ii]) # r32
  r[upper.tri(r)] <- 1-c(obj212$prob[ii], # r12
                         obj213$prob[ii], # r13
                         obj223$prob[ii]) # r23
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
  ### (2) Repeat (t = 1, 2, 3, ..., K, 1, ...)
  t <- 1
  while(TRUE){
    a <- 1/Q[t,t]
    b <- t(p) %*% Q %*% p
    p[t,] <- a * (as.vector(Q[t,]) %*% p  + b)
    
    ## normalize
    p <- p/sum(p)
    
    ## Condition (21) check
    tmp <- Q %*% p
    tmp2 <- outer(tmp, tmp, "-")
    tmp3 <- tmp2[upper.tri(tmp2)]
    if (length(unique(sign(tmp))) == 1 && abs(tmp3) < 1e-05) break
    
    ## re-indexing t
    if (t == K) {
      t <- 1
    } else {
      t <- (t + 1)
    }
  }
  
  # Save results
  result[[ii]] <-  p
}

# install.packages("doParallel")
library(doParallel)
# install.packages("foreach")
library(foreach)

# detectCores(all.tests = FALSE, logical = TRUE)

cl<-makeCluster(20)  # choose the number of cores
registerDoParallel(cl) # register clusters

system.time(
  mm<-foreach(ii=1:n.test) %dopar%
    working1(ii)
)

result.list<-lapply(1:12,function(j) 
  t(sapply(1:n.sim, function(k) mm[[k]][j,])))

result.list

# stopCluster(cl)

a <- as.vector(rep(0,30))
for (i in 1:30){
  a[i] <- which.max(result[[i]]$par)
}

a
test.y




  
  