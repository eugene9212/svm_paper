######################
## For class 4 Data ##
######################
rm(list = ls())

#### load packages & R code ####
library(mvtnorm)
library(fda)

setwd('C:/Users/eugene/Desktop/SVM/shared/R code/')
source('eu/data_gen/multiclass/class.4.R') # class4 generator

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
K <- 4
t <- seq(0, 1, by = 0.03)
L <- 10
error <- 0.1

n.train <- 120
n.test <- 40
n <- n.train + n.test

# Data generation
data <- class.4(n, error, t = t, seed = 1)

id <- sample(1:n, n.train)

train.x <- data$x[id]
train.y <- data$y[id]
test.x <- data$x[-id]
test.y <- data$y[-id]

# Divide train set pairwise
# Extract index
trn.idx12 <- which(train.y==1 | train.y==2)
trn.idx13 <- which(train.y==1 | train.y==3)
trn.idx14 <- which(train.y==1 | train.y==4)
trn.idx23 <- which(train.y==2 | train.y==3)
trn.idx24 <- which(train.y==2 | train.y==4)
trn.idx34 <- which(train.y==3 | train.y==4)
# train.x
train.x.12 <- train.x[trn.idx12]
train.x.13 <- train.x[trn.idx13]
train.x.14 <- train.x[trn.idx14]
train.x.23 <- train.x[trn.idx23]
train.x.24 <- train.x[trn.idx24]
train.x.34 <- train.x[trn.idx34]
# train.y
train.y.12 <- train.y[trn.idx12]
train.y.13 <- train.y[trn.idx13]
train.y.14 <- train.y[trn.idx14]
train.y.23 <- train.y[trn.idx23]
train.y.24 <- train.y[trn.idx24]
train.y.34 <- train.y[trn.idx34]

# Change y value into -1 and 1
train.y.12 <- ifelse(train.y.12 == 2, 1, -1)
train.y.13 <- ifelse(train.y.13 == 3, 1, -1)
train.y.14 <- ifelse(train.y.14 == 4, 1, -1)
train.y.23 <- ifelse(train.y.23 == 3, 1, -1)
train.y.24 <- ifelse(train.y.24 == 4, 1, -1)
train.y.34 <- ifelse(train.y.34 == 4, 1, -1)

####=======================     train data      ===============================####
# calculate the pi path
obj12 <- fsvm.prob(train.x.12, train.y.12, t, L)
obj13 <- fsvm.prob(train.x.13, train.y.13, t, L)
obj14 <- fsvm.prob(train.x.14, train.y.14, t, L)
obj23 <- fsvm.prob(train.x.23, train.y.23, t, L)
obj24 <- fsvm.prob(train.x.24, train.y.24, t, L)
obj34 <- fsvm.prob(train.x.34, train.y.34, t, L)

####======================= predict probability ===============================####
obj212 <- predict.fsvm.prob(obj12, test.x)
obj213 <- predict.fsvm.prob(obj13, test.x)
obj214 <- predict.fsvm.prob(obj14, test.x)
obj223 <- predict.fsvm.prob(obj23, test.x)
obj224 <- predict.fsvm.prob(obj24, test.x)
obj234 <- predict.fsvm.prob(obj34, test.x)

# # Save Data
# save(obj212, obj213, obj223, file="C:/Users/eugene/Desktop/class3_0.1.RData")

####============================= Pairwise Calculation =====================================####
###=============================== Algorithm 2. ============================================####
# # load Data
# load(file="class3_0.1.RData") # when error is 0.1
# load(file="class3_0.5.RData") # when error is 0.5
# load(file="class3_0.9.RData") # when error is 0.9
# K <- 3
# n <- 120
# n.test <- 30
# 
result <- as.list(1:n.test)

# Pairwise Coupling
for(ii in 1:n.test){
  
  # Create r matrix
  r <- matrix(0, K, K)
  r[lower.tri(r)] <- c(obj212$prob[ii], # r21
                       obj213$prob[ii], # r31
                       obj214$prob[ii], # r41
                       obj223$prob[ii], # r32
                       obj224$prob[ii], # r42
                       obj234$prob[ii]) # r43
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
  
  ### (2) Repeat (t = 1, 2, 3, ..., K, 1, ...)
  t <- 1
  iter.n <- 1
  
  while(TRUE){
    a <- 1/Q[t,t]
    b <- t(p) %*% Q %*% p
    p[t,] <- a * ( -as.vector(Q[t,-t]) %*% p[-t]  + b)
    
    ## normalize
    p <- p/sum(p)
    
    ## Condition (21) check
    tmp <- Q %*% p
    tmp2 <- matrix(outer(tmp, tmp, "-"), K, K)
    tmp3 <- tmp2[upper.tri(tmp2)]
    print(c(tmp3,iter.n))
    print(p)
    idx <- which(max(p) == p)
    c <- c(p[idx] - p[-idx])
    # if (length(unique(sign(tmp))) == 1 && abs(tmp3) < 1e-03)
    #   if(sum(p) == 1 && c > 0.99) break
    if (length(unique(sign(tmp))) == 1 && abs(tmp3) < 1e-05 && sum(p) == 1) break
    
    ## re-indexing t
    if (t == K) {
      t <- 1
    } else {
      t <- (t + 1)
    }
    
    iter.n <- iter.n + 1 # counting the iteration
  }
  
  # Save results
  result[[ii]] <- list(p = p, iter = iter.n, Q = Q)
}

result[[20]]$p
test.y

p.class <- c()
pred.p <- c()

for (ii in 1:n.test){
  p.class[ii] <- which(max(result[[ii]]$p) == result[[ii]]$p)
  pred.p[ii] <- max(result[[ii]]$p)
}

test.y
p.class
round(pred.p,2)



