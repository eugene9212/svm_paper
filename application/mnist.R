##########################################################################
# Goal : to compare mnist data by Bayes classifier, NN, flogistic, FSVM ##
##########################################################################
# install.packages("keras")
library(keras)
# install_keras()

###################
## Load the Data ##
###################
mnist_lst <- dataset_mnist()
str(mnist_lst) # Check the structure of data(MNIST)

## Load train_x,y
train_x <- mnist_lst$train$x
train_y <- mnist_lst$train$y

## Load test_x,y
test_x <- mnist_lst$test$x
test_y <- mnist_lst$test$y

## Data Flattening (2D array --> 1D array)
train_x <- array(as.numeric(train_x), dim = c(dim(train_x)[[1]], 784))
test_x <- array(as.numeric(test_x), dim = c(dim(test_x)[[1]], 784))

## AT this code, onlye digit 0 and 1 data are used
## Because of testing binary classification

## for 0,1 digit (train)
trn.ind0 <- which(train_y == 0) # extract index for digit0
trn.ind1 <- which(train_y == 1) # extract index for digit1
train_x0 <- train_x[trn.ind0,]  # create train.x for digit0 (5923 X 784)
train_y0 <- train_y[trn.ind0]
train_x1 <- train_x[trn.ind1,]  # create train.x for digit1 (6742 X 784)
train_y1 <- train_y[trn.ind1]

## for 0,1 digit (test)
tst.ind0 <- which(test_y == 0) # extract index for digit0
tst.ind1 <- which(test_y == 1) # extract index for digit1
test_x0 <- test_x[tst.ind0,]   # create test.x for digit0 (980 X 784)
test_y0 <- test_y[tst.ind0]
test_x1 <- test_x[tst.ind1,]   # create test.x for digit0 (1135 X 784) 
test_y1 <- test_y[tst.ind1]

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

##############################
## METHOD[1] Functional SVM ##
##############################
## Setting parameters
t <- seq(0, 1, by = 1/783)
L <- 10

## create train data
train_y0 <- ifelse(train_y0 == 0, -1, NA)
train.x <- rbind(train_x0, train_x1)
train.y <- rbind(as.matrix(train_y0,ncol=1), as.matrix(train_y1,ncol=1))
train.y <- as.vector(train.y) # length(train.y) = 12665

## Sampling train data(500) ##
num1 <- 500

## index
set.seed(1)
index <- sample(length(train.y), num1)
mean(train.y[index]) # whether it is balanced data

## sample of train.y
train.y <- train.y[index]

## sample of train.x
train.x <- train.x[index,]

u.list <- as.list(1:num1)
for (i in 1:dim(train.x)[1]){
  u.list[[i]] <- matrix(train.x[i,], nrow=1,byrow=TRUE)
}
trn.x.list <- u.list

####=======================     train data      ===============================####
obj <- fsvm.prob(trn.x.list, train.y, t, L)
####=======================     Model training completed    ===================####

## create test data 
test_y0 <- ifelse(test_y0 == 0, -1, NA)
test.x <- rbind(test_x0, test_x1)
test.y <- rbind(as.matrix(test_y0,ncol=1), as.matrix(test_y1,ncol=1))
test.y <- as.vector(test.y) # length(test.y) = 2115

## Sampling test data (200)
num2 <- 200

## index
set.seed(3)
index <- sample(length(test.y), num2)
mean(test.y[index]) # whether it is balanced data

## sample of test.y
test.y <- test.y[index]

## sample of test.x
test.x <- test.x[index,]

u.list <- as.list(1:num2)
for (i in 1:dim(test.x)[1]){
  u.list[[i]] <- matrix(test.x[i,], nrow=1,byrow=TRUE)
}
tst.x.list <- u.list

####======================= predict probability ===============================####
obj2 <- predict.fsvm.prob(obj, tst.x.list)
prob <- obj2$prob

# Boxplot of pi.star 
boxplot(prob[test.y == 1], prob[test.y == -1], xlab=paste("SVM"), ylim=c(0,1))

# Criteria
CRE <- -1/length(test.y)*(sum(1/2*(1+test.y)*log(prob) + 1/2*(1-test.y)*log(1-prob)))
CRE

###################################
## METHOD[2] Functional Logistic ##
###################################
## load packages
library(mvtnorm); library(fda); library(fda.usc)

## Setting parameters
t <- seq(0, 1, by = 1/783)
rangeval <- quantile(t, c(0,1))

#### Transform train.x and train.y
X <- trn.x.list
x.matrix <- matrix(unlist(X), nrow = length(X), byrow = T)
D <- list("data" = x.matrix, "argvals" = t, rangeval = rangeval)
attr(D, "class") <- "fdata"

Y <- train.y
index <- Y < 0
Y[index] <- 0
dataf <- as.data.frame(Y)

basis.obj <- create.bspline.basis(rangeval, 10) # L : number of basis
basis.x <- list("x"=basis.obj)
f <- Y ~ x
ldata <- list("df"=dataf,"x"=D)
res <- fregre.glm(f, family=binomial(link = "logit"), data=ldata, basis.x=basis.x)

# Predict
test.X <- tst.x.list
test.x.matrix <- matrix(unlist(test.X), nrow = length(test.X), byrow = T)
test.D <- list("data" = test.x.matrix, "argvals" = t, rangeval = rangeval)
attr(test.D, "class") <- "fdata"

test.Y <- test.y
index <- test.Y < 0
test.Y[index] <- 0
test.dataf <- as.data.frame(test.Y)

f <- test.Y ~ x
newldata <- list("df"=test.dataf,"x"=test.D)

# predict(res, newldata, type = "response")
pred.glm <- predict.fregre.glm(res, newldata)
prob <- pred.glm

boxplot(prob[test.Y == 1], prob[test.Y == 0], xlab=paste("logit"), ylim=c(0,1))

# Criteria
CRE <- -1/length(test.Y)*(sum(test.Y*log(prob) + (1-test.Y)*log(1-prob)))
CRE




