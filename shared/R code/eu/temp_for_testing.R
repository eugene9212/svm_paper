rm(list = ls())

#### load packages & R code ####
library(mvtnorm)
library(fda)

setwd('C:/Users/eugene/Desktop/SVM/shared/R code/')
source('eu/data_gen/multiclass/linear.par.K.R')    # linear with parallel        (beta)
source('eu/data_gen/multiclass/gp.1dim.I.K.R')

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
t <- seq(0, 1, by = 0.05)
L <- 10
beta <- 1
K <- 3

n.train <- 90
n.test <- 30
n <- n.train + n.test

# Data generation
# data <- linear.par3(n, beta, t = t, seed = 1)
data <- gp.1dim.I.K(n, beta, K, t = t, seed = 1)
# plot data_gen
idx1 <- which(data$y == 1)
idx2 <- which(data$y == 2)
idx3 <- which(data$y == 3)

a <- c()
for (i in 1:n) a <- c(a, unlist(data$x[[i]]))
plot(t, data$x[[1]], ylim = c(min(a),max(a)))
for (i in idx1[1:20]) lines(t, data$x[[i]], col = 1)
for (i in idx2[1:20]) lines(t, data$x[[i]], col = 2)
for (i in idx3[1:20]) lines(t, data$x[[i]], col = 3)

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
  
  result[[ii]] <-  p
}

a <- as.vector(rep(0,30))
for (i in 1:30){
  a[i] <- which.max(result[[i]]$par)
}

a
test.y





  
  