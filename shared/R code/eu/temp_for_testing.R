rm(list = ls())

#### load packages & R code ####
library(mvtnorm)
library(fda)

setwd('C:/Users/eugene/Desktop/SVM/shared/R code/')
source('eu/data_gen/multiclass/linear.par3.R')    # linear with parallel        (beta)
source('eu/data_gen/multiclass/gp.1dim.I3.R')

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

n.train <- 100
n.test <- 30
n <- n.train + n.test

# Data generation
# data <- linear.par3(n, beta, t = t, seed = 1)
data <- gp.1dim.I3(n, beta, t = t, seed = 1)
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

# Divide test set pairwise
# tst.idx12 <- which(test.y==1 | test.y==2)
# tst.idx13 <- which(test.y==1 | test.y==3)
# tst.idx23 <- which(test.y==2 | test.y==3)
# test.x.12 <- test.x[tst.idx12]
# test.x.13 <- test.x[tst.idx13]
# test.x.23 <- test.x[tst.idx23]
# test.y.12 <- test.y[tst.idx12]
# test.y.13 <- test.y[tst.idx13]
# test.y.23 <- test.y[tst.idx23]

####=======================     train data      ===============================####
# calculate the pi path
obj12 <- fsvm.prob(train.x.12, train.y.12, t, L)
obj13 <- fsvm.prob(train.x.13, train.y.13, t, L)
obj23 <- fsvm.prob(train.x.23, train.y.23, t, L)

####======================= predict probability ===============================####
obj212 <- predict.fsvm.prob(obj12, test.x)
obj213 <- predict.fsvm.prob(obj13, test.x)
obj223 <- predict.fsvm.prob(obj23, test.x)

result <- as.list(1:30)
for (i in 1:30){
  if (obj212$prob[i] < 0.5){
    p1 <- 1 - obj212$prob[i]
    p2 <- 1 - p1
  } else {
    p2 <- obj212$prob[i]
    p1 <- 1 - p2
  }
  
  r12 <- p1/(p1+p2) ; r21 <- c(1 - r12)
  
  if (obj13$prob[i] < 0.5){
    p1 <- 1 - obj213$prob[i]
    p3 <- 1 - p1
  } else {
    p3 <- obj213$prob[i]
    p1 <- 1 - p3
  }
  
  r13 <- p1/(p1+p3) ; r31 <- c(1 - r13)
  
  if (obj223$prob[i] < 0.5){
    p2 <- 1 - obj223$prob[i]
    p3 <- 1 - p2
  } else {
    p3 <- obj223$prob[i]
    p2 <- 1 - p3
  }
 
  r23 <- p2/(p2+p3) ; r32 <- c(1 - r23)

  # Pairwise Calculation
  opt_func <- function(p) {
    (r21*p[1] - r12*p[2])^2 + (r31*p[1] - r13*p[3])^2 +
      (r12*p[2] - r21*p[1])^2 + (r32*p[2] - r23*p[3])^2 +
      (r13*p[3] - r31*p[1])^2 + (r23*p[3] - r32*p[2])^2
  }
  result[[i]] <-  constrOptim(c(0.4,0.3,0.3), opt_func, ui=rbind(c(1, 1, 1), -c(1, 1, 1)), 
                              ci=c(1-1e-4,-1-1e-4), method='Nelder-Mead')
}

a <- as.vector(rep(0,30))
for (i in 1:30){
  a[i] <- which.max(result[[i]]$par)
}

a
test.y





# library(Rsolnp)
# opt_func <- function(p) {
#   (r21*p[1] + r31*p[1] - r12*p[2] - r13*p[3])^2 + 
#   (r12*p[2] + r32*p[2] - r21*p[1] - r23*p[3])^2 + 
#   (r13*p[3] + r23*p[3] - r31*p[1] - r32*p[2])^2
# }
# #specify the equality function. The number 15 (to which the function is equal)
# #is specified as an additional argument
# equal <- function(p) {
#   p[1] + p[2] + p[3] 
# }
# 
# #the optimiser - minimises by default
# solnp(c(0.4,0.3,0.3), #starting values (random - obviously need to be positive and sum to 15)
#       opt_func, #function to optimise
#       eqfun=equal, #equality function 
#       eqB=1,   #the equality constraint
#       LB=c(0,0,0), #lower bound for parameters i.e. greater than zero
#       # UB=c(100,100,100) #upper bound for parameters (I just chose 100 randomly)
#       )





  
  