rm(list = ls())
# TRial #1. Fourier basis

library(mvtnorm)
library(fda)
setwd('C:/Users/eugene/Desktop/SVM_R/shared/R code/')

source('fn/model1.R')
source('fn/multi.dim.fsvm.R')
source('fn/predict.multi.dim.fsvm.R')
source('fn/predictmultidimfsvm_fourier.R')
source('fn/fsvm.pi.path.R')
source('fn/fsvm.sub.pi.path.R')
source('fn/fsvm.lambda.path.R')
# trial
source('fn/multidimfsvm_fourier.R')

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir('KernSurf/R')


n <- 20      # of training sets

# generating data
beta <- 0.2 # controls distance between + and -
t <- seq(0, 1, by = 0.05) # time grid
n.r <- 4 # of rows in the lattice
n.c <- 4 # of columns in the lattice

seed <- 3 # seed for random numbrer

data <- model1(n, beta, t, n.r, n.c, seed)

x <- data$x
y <- data$y
t <- data$t
spatial.index <- data$spatial.index

n.grid <- nrow(spatial.index)
n.t <- length(t)

# Fitting multi dimensional functional svm
lambda <- .1 # regulrization parameter
rho <- 0.5 # tuning param for spatial correlation 
L <- 9  # number of Fourier basis
max.dist <- sqrt(2) + 1.0e-8 # maximum distance that defines neighborhood
period <- 2

dyn.load("KernSurf/temp/wsvmqp.dll")
obj <- multi.dim.fsvm_fourier(y, x, spatial.index, L, period, lambda, rho, max.dist) # 이게 제대로 되는건가
K<- obj$K # K : 20*20

# pi path가 잘 이해안됌
obj_pi <- fsvm.pi.path(lambda = 2, x, y, K)
str(obj_pi)

pi <- obj_pi$pi # 40
alpha <- obj_pi$alpha # 20 * 40 matrix
round(alpha,2)

# pi path
plot(0,0, type = "l", xlim = c(0,1), ylim = c(0,1),
     xlab = "pi", ylab = "alpha")
abline(v = pi, lty = 2, col = "gray")
for (ii in 1:n) lines(pi, alpha[ii,], col = 2)


#============================ <test> =============================================

## Test Set ##
# make prediction
test.n <- 20 # of test sets

# generating data
beta <- 0.2 # controls distance between + and -
t <- seq(0, 1, by = 0.05) # time grid
n.r <- 4 # of rows in the lattice
n.c <- 4 # of columns in the lattice

seed <- 3 # seed for random numbrer

test.data <- model1(test.n, beta, t, n.r, n.c, seed)

test.x <- test.data$x
test.y <- test.data$y
t <- test.data$t
spatial.index <- test.data$spatial.index

n.grid <- nrow(spatial.index)
n.t <- length(t)

obj1 <- predict.multi.dim.fsvm_fourier(obj,test.x) # obj는 train 에서 나온거.
obj1$new.fx # 20 outputs

# pi path
new.K <- obj1$new.K # K : 20*20 / train과 test의 K (근데 어디가 앞이고 뒤인지 모르겠어)
obj1_pi <- fsvm.pi.path(lambda = 2, test.x, test.y, new.K) # y는 train해야? test해야?

pi <- obj1_pi$pi # 39
alpha <- obj1_pi$alpha # 20 * 39 matrix
alpha0 <- obj1_pi$alpha0 # 39
round(alpha,2)

new.gx <- (new.K %*% (alpha * y)) # 20*39
new.fx <- matrix(0, test.n, length(pi))
pi.star <- rep(0, test.n)

# pi.star
for (ii in 1:test.n) {
  new.fx[ii,1] <- (alpha0[1] + new.gx[ii,1])/lambda
  
  for (jj in 2:length(pi)){
    
    new.fx[ii,jj] <- (alpha0[jj] + new.gx[ii,jj])/lambda
    a <- sign(new.fx[ii,jj-1])
    b <- sign(new.fx[ii,jj])
    
    if (a != b){
      pi.star[ii] <- (pi[jj-1] + pi[jj])/2
    }
  }
}
pi.star

# Boxplot of pi.star 
boxplot(pi.star[test.y == 1], pi.star[test.y != 1])

