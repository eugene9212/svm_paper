rm(list = ls())

#### load packages & R code ####
library(mvtnorm); library(fda); library(fda.usc)
setwd('C:/Users/eugene/Desktop/SVM/shared/R code/')
source('eu/data_gen/gp.1dim.I.R')     # GP with identity covariance (beta)
source('eu/data_gen/gp.1dim.cov1.R')  # GP with trivial covariance  (beta)
source('eu/data_gen/gp.1dim.sc.R')    # GP with trivial covariance  (beta)
source('eu/data_gen/gp.1dim.AR.R')    # GP with trivial covariance  (beta)
source('eu/data_gen/sc.R')            # non-linear with cross       (no beta)   
source('eu/data_gen/ss.R')            # non-linear with parallel    (beta)
source('eu/data_gen/linear.cross.R')  # linear with cross           (no beta)
source('eu/data_gen/linear.par.R')    # linear with parallel        (beta)

# set up
n.sim <- 50
t <- seq(0, 1, by = 0.05)
rangeval <- quantile(t, c(0,1))
L <- 10
beta <- 0.5

# storage
CRE.result<-matrix(0, n.sim, 1)
pi.result <- as.list(1:n.sim)

####========================= Simluation ==================================####
for (iter in 1:n.sim) {
  # iter<-37
  n.train <- 50
  n.test <- 40
  n <- n.train + n.test
  
  # Data generation (6 methods)
  set.seed(iter)
  # data <- gp.1dim.I(n, beta, t = t, seed = iter)
  # data <- gp.1dim.cov1(n, beta, t = t, seed = iter)
  # data <- gp.1dim.sc(n, t = t, seed = iter)
  # data <- gp.1dim.AR(n, beta, p = 5, rho = 0.5, t = t, seed = iter)
  # data <- sc(n, t = t, seed = iter)
  # data <- ss(n, beta, t = t, seed = iter)
  data <- linear.cross(n, t = t, seed = iter)
  # data <- linear.par(n, beta, t = t, seed = iter)
  
  id <- sample(1:n, n.train)
  
  train.x <- data$x[id]
  train.y <- data$y[id]
  test.x <- data$x[-id]
  test.y <- data$y[-id]
  
  print(iter)
  
  #### Transform train.x and train.y
  X <- train.x
  x.matrix <- matrix(unlist(X), nrow = length(X), byrow = T)
  D <- list("data" = x.matrix, "argvals" = t, rangeval = rangeval)
  attr(D, "class") <- "fdata"
  
  Y <- train.y; index <- Y < 0; Y[index] = 0
  dataf <- as.data.frame(Y)
  
  basis.obj <- create.bspline.basis(rangeval, 10) # L : number of basis
  basis.x <- list("x"=basis.obj)
  f <- Y ~ x
  ldata <- list("df"=dataf,"x"=D)
  res <- fregre.glm(f, family=binomial(link = "logit"), data=ldata, basis.x=basis.x)
  
  # Predict
  test.X <- test.x
  test.x.matrix <- matrix(unlist(test.X), nrow = length(test.X), byrow = T)
  test.D <- list("data" = test.x.matrix, "argvals" = t, rangeval = rangeval)
  attr(test.D, "class") <- "fdata"
  
  test.Y <- test.y; index <- test.Y<0; test.Y[index] = 0
  test.dataf <- as.data.frame(test.Y)
  
  f <- test.Y ~ x
  newldata <- list("df"=test.dataf,"x"=test.D)
  
  # predict(res, newldata, type = "response")
  pred.glm <- predict.fregre.glm(res, newldata)
  prob <- pred.glm
  
  boxplot(prob[test.Y == 1], prob[test.Y == 0], xlab=paste("logit",iter), ylim=c(0,1))
  
  # Criteria
  CRE <- -1/length(test.Y)*(sum(test.Y*log(prob) + (1-test.Y)*log(1-prob)))
  
  pi.result[[iter]] <- prob
  
  # results  
  CRE.result[iter,] <- c(CRE) 
}

mean(CRE.result)
# median(CRE.result)
# sum(CRE.result)



