# Goal : to compare mnist data by Bayes classifier, NN, flogistic, FSVM
# install.packages("keras")
library(keras)
# install_keras()
mnist_lst <- dataset_mnist()

str(mnist_lst)

# train_x,y
train_x <- mnist_lst$train$x
train_y <- mnist_lst$train$y

# test_x,y
test_x <- mnist_lst$test$x
test_y <- mnist_lst$test$y

# 4. 데이터 전처리 ----------------------------------
## 데이터 정규화(2D 배열 --> 1D 배열)
train_x <- array(as.numeric(train_x), dim = c(dim(train_x)[[1]], 784))
test_x <- array(as.numeric(test_x), dim = c(dim(test_x)[[1]], 784))

## for 0,1 digit (train)
trn.ind0 <- which(train_y == 0) ; trn.ind1 <- which(train_y == 1)

train_x0 <- train_x[trn.ind0,] ; train_y0 <- train_y[trn.ind0]
dim(train_x0) # 5923 X 784
train_x1 <- train_x[trn.ind1,] ; train_y1 <- train_y[trn.ind1]
dim(train_x1) # 6742 X 784

tst.ind0 <- which(test_y == 0) ; tst.ind1 <- which(test_y == 1)
test_x0 <- test_x[tst.ind0,] ; test_y0 <- test_y[tst.ind0]
dim(test_x0) # 980 X 784
test_x1 <- test_x[tst.ind1,] ; test_y1 <- test_y[tst.ind1]
dim(test_x1) # 1135 X 784

## plot as functional data(train)
# for digit=0
plot(train_x0[1,])
for (i in 1:dim(train_x0)[1]/2) lines(train_x0[i,], col = 1)

# for digit=1
for (i in 1:dim(train_x1)[1]/2) lines(train_x1[i,], col = 4)

## plot as functional data(test)
# for digit=0
plot(test_x0[1,])
for (i in 1:dim(test_x0)[1]/2) lines(test_x0[i,], col = 1)

# for digit=1
for (i in 1:dim(test_x1)[1]/2) lines(test_x1[i,], col = 4)

######################
## load the methods ##
######################
# method1. Functional SVM
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

t <- seq(0, 1, by = 1/783); L <- 10
# create train data
dim(train_x0)
train_y0 <- ifelse(train_y0 == 0, -1, NA)
train.x <- rbind(train_x0, train_x1)
train.y <- rbind(as.matrix(train_y0,ncol=1), as.matrix(train_y1,ncol=1))
train.y <- as.vector(train.y)
dim(train.y)
class(train.y)

# u.list <- as.list(1:dim(train.x)[1])
# for (i in 1:dim(train.x)[1]){
#   u.list[[i]] <- matrix(train.x[i,], nrow=1,byrow=TRUE)
# }

num <- 500
as.list <- as.list(1:num)
index <- sample(length(train.x.list),500)

for (i in index){
  u.list[[i]] <- matrix(train.x[i,], nrow=1,byrow=TRUE)
}
u.list[[1]]
train.y <- train.y[index]
train.y <- as.vector(train.y)
dim(train.y)
class(train.y)
str(train.y)
outer(train.y,train.y)
####=======================     train data      ===============================####
obj <- fsvm.prob(u.list, train.y, t, L)
######################==============ERROR==================#########################

####======================= predict probability ===============================####
obj2 <- predict.fsvm.prob(obj, test.x)
prob <- obj2$prob

# Boxplot of pi.star 
boxplot(prob[test.y == 1], prob[test.y == -1], xlab=paste("SVM",iter), ylim=c(0,1))

# Criteria
# delta = 1.0e-8
CRE <- -1/length(test.y)*(sum(1/2*(1+test.y)*log(prob) + 1/2*(1-test.y)*log(1-prob)))

