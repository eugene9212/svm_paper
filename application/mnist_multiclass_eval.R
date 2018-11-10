# Testing multiclass probability estimation #
# Data : MNIST #
rm(list = ls())

library(keras)
##########
## Data ##
##########
# Load the Data 
mnist_lst <- dataset_mnist()

## Load train_x,y
train_x <- mnist_lst$train$x
train_y <- mnist_lst$train$y

## Load test_x,y
test_x <- mnist_lst$test$x
test_y <- mnist_lst$test$y

## Data Flattening (2D array --> 1D array)
train_x <- 3*array(as.numeric(train_x), dim = c(dim(train_x)[[1]], 784))
test_x <- 3*array(as.numeric(test_x), dim = c(dim(test_x)[[1]], 784))

# Sampling train set
set.seed(1)
index <- sample(dim(train_x)[1], size = 1000, replace=FALSE)
train.x <- train_x[index,]
train.y <- train_y[index]

# Sampling test set
set.seed(1)
index1 <- sample(dim(test_x)[1], size = 500, replace=FALSE)
test.x <- test_x[index1,]
test.y <- test_y[index1]


# Divide train set pairwise --------------------------------------------------------------------(TRAIN)
# (1) Extract index
K <- 9
for (i in 0:(K-1)){
  for (j in (i+1):K){
    eval(parse(text=paste0("trn.idx",i,j," <- which(train.y==",i," | train.y==",j,")")))
    cat(i,j,' pair is clear\n')
  }
}

# (2) Split train.x pairwise

for (i in 0:(K-1)){
  for (j in (i+1):K){
    eval(parse(text=paste0("train.x",i,j," <- train.x[trn.idx",i,j,",]")))
    cat(i,j,' pair is clear\n')
  }
}


# (3) Split train.y pairwise
for (i in 0:(K-1)){
  for (j in (i+1):K){
    eval(parse(text=paste0("train.y",i,j," <- train.y[trn.idx",i,j,",]")))
    cat(i,j,' pair is clear\n')
  }
}

# (4) Change y value into -1 and 1
for (i in 0:(K-1)){
  for (j in (i+1):K){
    eval(parse(text=paste0("train.y",i,j," <- ifelse(train.y",i,j," == ",j,", 1, -1)")))
    cat(i,j,' pair is clear\n')
  }
}
# Divide test set pairwise --------------------------------------------------------------------(TEST)
# (1) Extract index (test)
K <- 9
for (i in 0:(K-1)){
  for (j in (i+1):K){
    eval(parse(text=paste0("tst.idx",i,j," <- which(test.y==",i," | test.y==",j,")")))
    cat(i,j,' pair is clear\n')
  }
}

# (2) Split test.x pairwise

for (i in 0:(K-1)){
  for (j in (i+1):K){
    eval(parse(text=paste0("test.x",i,j," <- test.x[tst.idx",i,j,",]")))
    cat(i,j,' pair is clear\n')
  }
}

# (3) Split train.y pairwise
for (i in 0:(K-1)){
  for (j in (i+1):K){
    eval(parse(text=paste0("test.y",i,j," <- test.y[tst.idx",i,j,",]")))
    cat(i,j,' pair is clear\n')
  }
}

# (4) Change y value into -1 and 1
for (i in 0:(K-1)){
  for (j in (i+1):K){
    eval(parse(text=paste0("test.y",i,j," <- ifelse(test.y",i,j," == ",j,", 1, -1)")))
    cat(i,j,' pair is clear\n')
  }
}

######################
## load the methods ##
######################
library(fda)
setwd('C:/Users/eugene/Desktop/SVM/shared/R code/')
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

# Parameters
L <- 10
t <- seq(0, 1, by = 1/783)

####=======================     train data      ===============================####
#### convert matrix into train.x.list ####
# class 0 pair
num <- dim(train.x01)[1]
train.x01.list <- as.list(1:num)
for (i in 1:num) train.x01.list[[i]] <- matrix(train.x01[i,], nrow=1,byrow=TRUE)

num <- dim(train.x02)[1]
train.x02.list <- as.list(1:num)
for (i in 1:num) train.x02.list[[i]] <- matrix(train.x02[i,], nrow=1,byrow=TRUE)

num <- dim(train.x03)[1]
train.x03.list <- as.list(1:num)
for (i in 1:num) train.x03.list[[i]] <- matrix(train.x03[i,], nrow=1,byrow=TRUE)

num <- dim(train.x04)[1]
train.x04.list <- as.list(1:num)
for (i in 1:num) train.x04.list[[i]] <- matrix(train.x04[i,], nrow=1,byrow=TRUE)

num <- dim(train.x05)[1]
train.x05.list <- as.list(1:num)
for (i in 1:num) train.x05.list[[i]] <- matrix(train.x05[i,], nrow=1,byrow=TRUE)

num <- dim(train.x06)[1]
train.x06.list <- as.list(1:num)
for (i in 1:num) train.x06.list[[i]] <- matrix(train.x06[i,], nrow=1,byrow=TRUE)

num <- dim(train.x07)[1]
train.x07.list <- as.list(1:num)
for (i in 1:num) train.x07.list[[i]] <- matrix(train.x07[i,], nrow=1,byrow=TRUE)

num <- dim(train.x08)[1]
train.x08.list <- as.list(1:num)
for (i in 1:num) train.x08.list[[i]] <- matrix(train.x08[i,], nrow=1,byrow=TRUE)

num <- dim(train.x09)[1]
train.x09.list <- as.list(1:num)
for (i in 1:num) train.x09.list[[i]] <- matrix(train.x09[i,], nrow=1,byrow=TRUE)

# class 1 pair
num <- dim(train.x12)[1]
train.x12.list <- as.list(1:num)
for (i in 1:num) train.x12.list[[i]] <- matrix(train.x12[i,], nrow=1,byrow=TRUE)

num <- dim(train.x13)[1]
train.x13.list <- as.list(1:num)
for (i in 1:num) train.x13.list[[i]] <- matrix(train.x13[i,], nrow=1,byrow=TRUE)

num <- dim(train.x14)[1]
train.x14.list <- as.list(1:num)
for (i in 1:num) train.x14.list[[i]] <- matrix(train.x14[i,], nrow=1,byrow=TRUE)

num <- dim(train.x15)[1]
train.x15.list <- as.list(1:num)
for (i in 1:num) train.x15.list[[i]] <- matrix(train.x15[i,], nrow=1,byrow=TRUE)

num <- dim(train.x16)[1]
train.x16.list <- as.list(1:num)
for (i in 1:num) train.x16.list[[i]] <- matrix(train.x16[i,], nrow=1,byrow=TRUE)

num <- dim(train.x17)[1]
train.x17.list <- as.list(1:num)
for (i in 1:num) train.x17.list[[i]] <- matrix(train.x17[i,], nrow=1,byrow=TRUE)

num <- dim(train.x18)[1]
train.x18.list <- as.list(1:num)
for (i in 1:num) train.x18.list[[i]] <- matrix(train.x18[i,], nrow=1,byrow=TRUE)

num <- dim(train.x19)[1]
train.x19.list <- as.list(1:num)
for (i in 1:num) train.x19.list[[i]] <- matrix(train.x19[i,], nrow=1,byrow=TRUE)

# class 2 pair
num <- dim(train.x23)[1]
train.x23.list <- as.list(1:num)
for (i in 1:num) train.x23.list[[i]] <- matrix(train.x23[i,], nrow=1,byrow=TRUE)
num <- dim(train.x24)[1]
train.x24.list <- as.list(1:num)
for (i in 1:num) train.x24.list[[i]] <- matrix(train.x24[i,], nrow=1,byrow=TRUE)
num <- dim(train.x25)[1]
train.x25.list <- as.list(1:num)
for (i in 1:num) train.x25.list[[i]] <- matrix(train.x25[i,], nrow=1,byrow=TRUE)
num <- dim(train.x26)[1]
train.x26.list <- as.list(1:num)
for (i in 1:num) train.x26.list[[i]] <- matrix(train.x26[i,], nrow=1,byrow=TRUE)
num <- dim(train.x27)[1]
train.x27.list <- as.list(1:num)
for (i in 1:num) train.x27.list[[i]] <- matrix(train.x27[i,], nrow=1,byrow=TRUE)
num <- dim(train.x28)[1]
train.x28.list <- as.list(1:num)
for (i in 1:num) train.x28.list[[i]] <- matrix(train.x28[i,], nrow=1,byrow=TRUE)
num <- dim(train.x29)[1]
train.x29.list <- as.list(1:num)
for (i in 1:num) train.x29.list[[i]] <- matrix(train.x29[i,], nrow=1,byrow=TRUE)
# class 3 pair
num <- dim(train.x34)[1]
train.x34.list <- as.list(1:num)
for (i in 1:num) train.x34.list[[i]] <- matrix(train.x34[i,], nrow=1,byrow=TRUE)
num <- dim(train.x35)[1]
train.x35.list <- as.list(1:num)
for (i in 1:num) train.x35.list[[i]] <- matrix(train.x35[i,], nrow=1,byrow=TRUE)
num <- dim(train.x36)[1]
train.x36.list <- as.list(1:num)
for (i in 1:num) train.x36.list[[i]] <- matrix(train.x36[i,], nrow=1,byrow=TRUE)
num <- dim(train.x37)[1]
train.x37.list <- as.list(1:num)
for (i in 1:num) train.x37.list[[i]] <- matrix(train.x37[i,], nrow=1,byrow=TRUE)
num <- dim(train.x38)[1]
train.x38.list <- as.list(1:num)
for (i in 1:num) train.x38.list[[i]] <- matrix(train.x38[i,], nrow=1,byrow=TRUE)
num <- dim(train.x39)[1]
train.x39.list <- as.list(1:num)
for (i in 1:num) train.x39.list[[i]] <- matrix(train.x39[i,], nrow=1,byrow=TRUE)
# class 4 pair
num <- dim(train.x45)[1]
train.x45.list <- as.list(1:num)
for (i in 1:num) train.x45.list[[i]] <- matrix(train.x45[i,], nrow=1,byrow=TRUE)
num <- dim(train.x46)[1]
train.x46.list <- as.list(1:num)
for (i in 1:num) train.x46.list[[i]] <- matrix(train.x46[i,], nrow=1,byrow=TRUE)
num <- dim(train.x47)[1]
train.x47.list <- as.list(1:num)
for (i in 1:num) train.x47.list[[i]] <- matrix(train.x47[i,], nrow=1,byrow=TRUE)
num <- dim(train.x48)[1]
train.x48.list <- as.list(1:num)
for (i in 1:num) train.x48.list[[i]] <- matrix(train.x48[i,], nrow=1,byrow=TRUE)
num <- dim(train.x49)[1]
train.x49.list <- as.list(1:num)
for (i in 1:num) train.x49.list[[i]] <- matrix(train.x49[i,], nrow=1,byrow=TRUE)
# class 5 pair
num <- dim(train.x56)[1]
train.x56.list <- as.list(1:num)
for (i in 1:num) train.x56.list[[i]] <- matrix(train.x56[i,], nrow=1,byrow=TRUE)
num <- dim(train.x57)[1]
train.x57.list <- as.list(1:num)
for (i in 1:num) train.x57.list[[i]] <- matrix(train.x57[i,], nrow=1,byrow=TRUE)
num <- dim(train.x58)[1]
train.x58.list <- as.list(1:num)
for (i in 1:num) train.x58.list[[i]] <- matrix(train.x58[i,], nrow=1,byrow=TRUE)
num <- dim(train.x59)[1]
train.x59.list <- as.list(1:num)
for (i in 1:num) train.x59.list[[i]] <- matrix(train.x59[i,], nrow=1,byrow=TRUE)
# class 6 pair
num <- dim(train.x67)[1]
train.x67.list <- as.list(1:num)
for (i in 1:num) train.x67.list[[i]] <- matrix(train.x67[i,], nrow=1,byrow=TRUE)
num <- dim(train.x68)[1]
train.x68.list <- as.list(1:num)
for (i in 1:num) train.x68.list[[i]] <- matrix(train.x68[i,], nrow=1,byrow=TRUE)
num <- dim(train.x69)[1]
train.x69.list <- as.list(1:num)
for (i in 1:num) train.x69.list[[i]] <- matrix(train.x69[i,], nrow=1,byrow=TRUE)
# class 7 pair
num <- dim(train.x78)[1]
train.x78.list <- as.list(1:num)
for (i in 1:num) train.x78.list[[i]] <- matrix(train.x78[i,], nrow=1,byrow=TRUE)
num <- dim(train.x79)[1]
train.x79.list <- as.list(1:num)
for (i in 1:num) train.x79.list[[i]] <- matrix(train.x79[i,], nrow=1,byrow=TRUE)
# class 8 pair
num <- dim(train.x89)[1]
train.x89.list <- as.list(1:num)
for (i in 1:num) train.x89.list[[i]] <- matrix(train.x89[i,], nrow=1,byrow=TRUE)

#### calculate the pi path ####
# class 0 pair
obj01 <- fsvm.prob(train.x01.list, c(train.y01), t, L)
obj02 <- fsvm.prob(train.x02.list, c(train.y02), t, L)
obj03 <- fsvm.prob(train.x03.list, c(train.y03), t, L)
obj04 <- fsvm.prob(train.x04.list, c(train.y04), t, L)
obj05 <- fsvm.prob(train.x05.list, c(train.y05), t, L)
obj06 <- fsvm.prob(train.x06.list, c(train.y06), t, L)
obj07 <- fsvm.prob(train.x07.list, c(train.y07), t, L)
obj08 <- fsvm.prob(train.x08.list, c(train.y08), t, L)
obj09 <- fsvm.prob(train.x09.list, c(train.y09), t, L)
# class 1 pair
obj12 <- fsvm.prob(train.x12.list, c(train.y12), t, L)
obj13 <- fsvm.prob(train.x13.list, c(train.y13), t, L) #########################
obj14 <- fsvm.prob(train.x14.list, c(train.y14), t, L)
obj15 <- fsvm.prob(train.x15.list, c(train.y15), t, L)
obj16 <- fsvm.prob(train.x16.list, c(train.y16), t, L)
obj17 <- fsvm.prob(train.x17.list, c(train.y17), t, L)
obj18 <- fsvm.prob(train.x18.list, c(train.y18), t, L)
obj19 <- fsvm.prob(train.x19.list, c(train.y19), t, L)
# class 2 pair
obj23 <- fsvm.prob(train.x23.list, c(train.y23), t, L)
obj24 <- fsvm.prob(train.x24.list, c(train.y24), t, L)
obj25 <- fsvm.prob(train.x25.list, c(train.y25), t, L)
obj26 <- fsvm.prob(train.x26.list, c(train.y26), t, L)
obj27 <- fsvm.prob(train.x27.list, c(train.y27), t, L)
obj28 <- fsvm.prob(train.x28.list, c(train.y28), t, L)
obj29 <- fsvm.prob(train.x29.list, c(train.y29), t, L)
# class 3 pair
obj34 <- fsvm.prob(train.x34.list, c(train.y34), t, L)
obj35 <- fsvm.prob(train.x35.list, c(train.y35), t, L)
obj36 <- fsvm.prob(train.x36.list, c(train.y36), t, L)
obj37 <- fsvm.prob(train.x37.list, c(train.y37), t, L)
obj38 <- fsvm.prob(train.x38.list, c(train.y38), t, L)
obj39 <- fsvm.prob(train.x39.list, c(train.y39), t, L)
# class 4 pair
obj45 <- fsvm.prob(train.x45.list, c(train.y45), t, L)
obj46 <- fsvm.prob(train.x46.list, c(train.y46), t, L)
obj47 <- fsvm.prob(train.x47.list, c(train.y47), t, L)
obj48 <- fsvm.prob(train.x48.list, c(train.y48), t, L)
obj49 <- fsvm.prob(train.x49.list, c(train.y49), t, L)
# class 5 pair
obj56 <- fsvm.prob(train.x56.list, c(train.y56), t, L)
obj57 <- fsvm.prob(train.x57.list, c(train.y57), t, L)
obj58 <- fsvm.prob(train.x58.list, c(train.y58), t, L)
obj59 <- fsvm.prob(train.x59.list, c(train.y59), t, L)
# class 6 pair
obj67 <- fsvm.prob(train.x67.list, c(train.y67), t, L)
obj68 <- fsvm.prob(train.x68.list, c(train.y68), t, L)
obj69 <- fsvm.prob(train.x69.list, c(train.y69), t, L)
# class 7 pair
obj78 <- fsvm.prob(train.x78.list, c(train.y78), t, L)
obj79 <- fsvm.prob(train.x79.list, c(train.y79), t, L)
# class 8 pair
obj89 <- fsvm.prob(train.x89.list, c(train.y89), t, L)

####======================= predict probability ===============================####
## convert matrix into list
# class 0 pair
num <- dim(test.x)[1]
test.x.list <- as.list(1:num)
for (i in 1:num) test.x.list[[i]] <- matrix(test.x[i,], nrow=1,byrow=TRUE)
# class 0 pair
obj201 <- predict.fsvm.prob(obj01, test.x.list)
obj202 <- predict.fsvm.prob(obj02, test.x.list)
obj203 <- predict.fsvm.prob(obj03, test.x.list)
obj204 <- predict.fsvm.prob(obj04, test.x.list)
obj205 <- predict.fsvm.prob(obj05, test.x.list)
obj206 <- predict.fsvm.prob(obj06, test.x.list)
obj207 <- predict.fsvm.prob(obj07, test.x.list)
obj208 <- predict.fsvm.prob(obj08, test.x.list)
obj209 <- predict.fsvm.prob(obj09, test.x.list)
# class 1 pair
obj212 <- predict.fsvm.prob(obj12, test.x.list)
obj213 <- predict.fsvm.prob(obj13, test.x.list)
obj214 <- predict.fsvm.prob(obj14, test.x.list)
obj215 <- predict.fsvm.prob(obj15, test.x.list)
obj216 <- predict.fsvm.prob(obj16, test.x.list)
obj217 <- predict.fsvm.prob(obj17, test.x.list)
obj218 <- predict.fsvm.prob(obj18, test.x.list)
obj219 <- predict.fsvm.prob(obj19, test.x.list)
# class 2 pair
obj223 <- predict.fsvm.prob(obj23, test.x.list)
obj224 <- predict.fsvm.prob(obj24, test.x.list)
obj225 <- predict.fsvm.prob(obj25, test.x.list)
obj226 <- predict.fsvm.prob(obj26, test.x.list)
obj227 <- predict.fsvm.prob(obj27, test.x.list)
obj228 <- predict.fsvm.prob(obj28, test.x.list)
obj229 <- predict.fsvm.prob(obj29, test.x.list)
# class 3 pair
obj234 <- predict.fsvm.prob(obj34, test.x.list)
obj235 <- predict.fsvm.prob(obj35, test.x.list)
obj236 <- predict.fsvm.prob(obj36, test.x.list)
obj237 <- predict.fsvm.prob(obj37, test.x.list)
obj238 <- predict.fsvm.prob(obj38, test.x.list)
obj239 <- predict.fsvm.prob(obj39, test.x.list)
# class 4 pair
obj245 <- predict.fsvm.prob(obj45, test.x.list)
obj246 <- predict.fsvm.prob(obj46, test.x.list)
obj247 <- predict.fsvm.prob(obj47, test.x.list)
obj248 <- predict.fsvm.prob(obj48, test.x.list)
obj249 <- predict.fsvm.prob(obj49, test.x.list)
# class 5 pair
obj256 <- predict.fsvm.prob(obj56, test.x.list)
obj257 <- predict.fsvm.prob(obj57, test.x.list)
obj258 <- predict.fsvm.prob(obj58, test.x.list)
obj259 <- predict.fsvm.prob(obj59, test.x.list)
# class 6 pair
obj267 <- predict.fsvm.prob(obj67, test.x.list)
obj268 <- predict.fsvm.prob(obj68, test.x.list)
obj269 <- predict.fsvm.prob(obj69, test.x.list)
# class 7 pair
obj278 <- predict.fsvm.prob(obj78, test.x.list)
obj279 <- predict.fsvm.prob(obj79, test.x.list)
# class 8 pair
obj289 <- predict.fsvm.prob(obj89, test.x.list)

####============================= Pairwise Calculation =====================================####
###=============================== Algorithm 2. ============================================####
n.test <- length(test.y)
result <- as.list(1:n.test)
K <- 10

# Pairwise Coupling
for(ii in 1:n.test){
  r <- matrix(0, K, K)
  r[lower.tri(r)] <- c(obj201$prob[ii], obj202$prob[ii], obj203$prob[ii], obj204$prob[ii], obj205$prob[ii], obj206$prob[ii], obj207$prob[ii], obj208$prob[ii], obj209$prob[ii], 
      obj212$prob[ii], obj213$prob[ii], obj214$prob[ii], obj215$prob[ii], obj216$prob[ii], obj217$prob[ii], obj218$prob[ii], obj219$prob[ii],
      obj223$prob[ii], obj224$prob[ii], obj225$prob[ii], obj226$prob[ii], obj227$prob[ii], obj228$prob[ii], obj229$prob[ii],
      obj234$prob[ii], obj235$prob[ii], obj236$prob[ii], obj237$prob[ii], obj238$prob[ii], obj239$prob[ii],
      obj245$prob[ii], obj246$prob[ii], obj247$prob[ii], obj248$prob[ii], obj249$prob[ii],
      obj256$prob[ii], obj257$prob[ii], obj258$prob[ii], obj259$prob[ii],
      obj267$prob[ii], obj268$prob[ii], obj269$prob[ii],
      obj278$prob[ii], obj279$prob[ii],
      obj289$prob[ii])
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
    # print(c(tmp3,iter.n))
    # print(p)
    idx <- which(max(p) == p)[1]
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

result[[1]]$p
test.y

p.class <- c()
pred.p <- c()

for (ii in 1:n.test){
  p.class[ii] <- which(max(result[[ii]]$p) == result[[ii]]$p)
  pred.p[ii] <- max(result[[ii]]$p)
}

dim(test.y)
length(p.class)
sum(test.y == p.class)
test.y[test.y == p.class]
pred.p





