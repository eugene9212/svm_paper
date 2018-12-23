# 18.12.23 #
#### To find out the misleading part
#### by changing the error level
#### and the number of train set

### Goal : To find out why an increase of n.train induces an increas of CRE criterion (which is bad)
rm(list = ls())

#### load packages & R code ####
library(mvtnorm)
library(fda)
library(fda.usc)

setwd('C:/Users/eugene/Desktop/SVM/shared/R code/')
source('eu/data_gen/binary/gp.1dim.sc.R')  # GP with trivial covariance  (beta)
source('eu/data_gen/binary/gp.1dim.ss.R')  # GP with trivial covariance  (beta)
source('eu/data_gen/binary/linear.cross.R')  # linear with cross           (no beta)
source('eu/data_gen/binary/linear.par.R')    # linear with parallel        (beta)

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
n.sim <- 30
t <- seq(0, 1, by = 0.05)
L <- 10
beta <- 1
rho <- 0.1

n.rho <- c(0.1,0.3,0.5,0.6)

# create a storage by method
result.by.gen <- as.list(1:4)
result <- matrix(0, nrow = length(n.train) * length(n.rho), ncol = 8)
colnames(result) <- c("svm.cre","fl.cre","svm.acc","fl.acc","diff.svm","diff.fl",
                      "w.diff.svm","w.diff.fl")
rownames(result) <- c("100/0.1","200/0.1","300/0.1","500/0.1",
                      "100/0.3","200/0.3","300/0.3","500/0.3",
                      "100/0.5","200/0.5","300/0.5","500/0.5",
                      "100/0.6","200/0.6","300/0.6","500/0.6")
# result


n.train <- c(100,200,300,500)
n.test <- 1000

cov <- "AR"


for(n.fn in 1:4){
  
  # error
  e <- c(0.5,1)
  # Storage by one method
  result <- matrix(0, nrow = length(n.train) * length(n.rho), ncol = 8)
  colnames(result) <- c("svm.cre","fl.cre","svm.acc","fl.acc","diff.svm","diff.fl","w.diff.svm","w.diff.fl")
  rownames(result) <- c("100/0.1","200/0.1","300/0.1","500/0.1",
                        "100/0.3","200/0.3","300/0.3","500/0.3",
                        "100/0.5","200/0.5","300/0.5","500/0.5",
                        "100/0.6","200/0.6","300/0.6","500/0.6")
  
  for(r in 1:length(n.rho)){
    rho <- n.rho[r]
    
    for(trn in 1:length(n.train)){
      n.trn <- n.train[trn]
      
      for(j in 2){
        error <- e[j]
        
        # storage
        ## for test.y
        ans <- as.list(1:n.sim)
        ans.p <- as.list(1:n.sim)
        ## for Functional logistic
        CRE.fl.result<-matrix(0, n.sim, 1)
        pi.fl.result <- as.list(1:n.sim)
        ## for SVM
        CRE.svm.result<-matrix(0, n.sim, 1)
        pi.svm.result <- as.list(1:n.sim)
        
        ####========================= Simluation ==================================####
        for (iter in 1:n.sim) {
          # iter<-1
          n <- n.trn + n.test
          
          # Data generation (4 methods with cov=Identity)
          set.seed(iter)
          if(n.fn == 1) {
            data <- gp.1dim.ss(n, error = error, beta, cov = cov, rho = rho, t = t, seed = iter)
          } else if (n.fn == 2) {
            data <- gp.1dim.sc(n, error = error, cov = cov, rho = rho, t = t, seed = iter)
          } else if (n.fn == 3) {
            data <- linear.cross(n, error = error, cov = cov, rho = rho, t = t, seed = iter)
          } else if (n.fn == 4) {
            data <- linear.par(n, error = error, beta, cov = cov, rho = rho, t = t, seed = iter)
          } else {
            warning("data generation fn is not matched")
          }
          
          id <- sample(1:n, n.trn)
          
          train.x <- data$x[id]
          train.y <- data$y[id]
          test.x <- data$x[-id]
          test.y <- data$y[-id]
          
          ans.p[[iter]] <- data$true.p[-id]
          
          print(iter)
          
          #========================================== Functional SVM ====================####
          ####=======================     train data
          svm.obj <- tryCatch({
            fsvm.prob(train.x, train.y, t, L)
          }, error = function(e){
            fsvm.prob(train.x, train.y, t, L, ridge=1e-3)
          })
          # Error in solve.default(Kstar, al)
          # svm.obj <- fsvm.prob(train.x, train.y, t, L, ridge) # calculate the pi path
          
          ####=======================     test data
          svm.obj2 <- predict.fsvm.prob(svm.obj, test.x)
          svm.prob <- svm.obj2$prob
          # cbind(test.y,svm.prob)
          
          #========================================== Functional logistic ===============####
          #### Transform train.x and train.y
          train.fy <- ifelse(train.y == 1, 1, 0)
          train.fy <- data.frame(train.fy)
          train.f.x.matrix <- matrix(unlist(train.x), nrow = length(train.x), byrow = T)
          train.fx <- fdata(train.f.x.matrix,argvals=t,rangeval=range(t))
          
          #### FDA
          nbasis.x=L; nbasis.b=L
          
          # create basis
          basis1=create.bspline.basis(rangeval=range(t),nbasis=nbasis.x)
          basis2=create.bspline.basis(rangeval=range(t),nbasis=nbasis.b)
          # Create basis n ldata before fitting the model
          basis.x=list("x"=basis1) ; basis.b=list("x"=basis2)
          
          # formula
          f=train.fy~x
          # as.factor
          train.fy <- train.fy$train.fy
          ldata=list("df"=train.fy,"x"=train.fx)
          
          ####=======================     train data
          # Fit the model
          fl.fit <- fregre.glm(f,familiy=binomial(link = "logit"), data=ldata, basis.x=basis.x, basis.b=basis.b, control =list(maxit=1000))
          
          ####=======================     test data
          # test.x -> test.fx(fdata)
          test.fx.matrix <- matrix(unlist(test.x), nrow = length(test.x), byrow = T)
          test.fx <- fdata(test.fx.matrix, argvals=t, rangeval=range(t))
          
          # create newldata
          newldata <- list("x"=test.fx)
          
          # predict
          prob <- predict.fregre.glm(fl.fit, newldata)
          prob <- ifelse(prob < 0,1e-10,prob)
          prob <- ifelse(prob > 1,1-1e-10,prob)
          
          # Predicted Probability of Funtional logistic
          pi.fl.result[[iter]] <- c(prob)
          # cbind(test.y, prob)
          
          # Criteria (CRE)
          delta <- 1e-10
          CRE.fl <- -1/length(test.y)*(sum(1/2*(1+test.y)*log(prob+delta) + 1/2*(1-test.y)*log(1-prob+delta)))
          CRE.fl.result[iter,] <- CRE.fl
          
          # Predicted Probability of FSVM
          pi.svm.result[[iter]] <- c(svm.prob)
          
          # Criteria (CRE)
          CRE.svm <- -1/length(test.y)*(sum(1/2*(1+test.y)*log(svm.prob+delta) + 1/2*(1-test.y)*log(1-svm.prob+delta)))
          CRE.svm.result[iter,] <- CRE.svm
          
          # Box Plot
          # boxplot(svm.prob[test.y == 1], svm.prob[test.y == -1], xlab=paste("svm",iter), ylim=c(0,1))
          # boxplot(prob[test.y == 1], prob[test.y == -1], xlab=paste("logit",iter), ylim=c(0,1))
          
          # store the answers
          ans[[iter]] <- test.y
        }
        # Critereon (1) Cross Entropy ======================================================
        fl.cre <- round(mean(CRE.fl.result),digits = 3)
        svm.cre <- round(mean(CRE.svm.result),digits = 3)
        
        # Critereon (2) Accuracy ===========================================================
        svm <- matrix(0, ncol = n.sim)
        flog <- matrix(0, ncol = n.sim)
        for(k in 1:n.sim){
          flog[,k] <- sum(ifelse(pi.fl.result[[k]]<=0.5, -1, 1) == ans[[k]])
          svm[,k] <- sum(ifelse(pi.svm.result[[k]]<=0.5, -1, 1) == ans[[k]])
        }
        svm.acc <- round(mean(svm/n.test), digits = 3)
        fl.acc <- round(mean(flog/n.test), digits = 3)
        
        # Critereon (3) Distance btw true p & hat p =========================================
        predict.p.svm <- pi.svm.result
        predict.p.fl <- pi.fl.result
        
        for(i in 1:n.sim){
          # i<-1
          class.fl <- c(); idx.fl <- c(); class.svm <- c(); idx.svm <- c()
          
          class.fl <- ifelse(pi.fl.result[[i]]<=0.5, -1, 1)
          idx.fl <- which(class.fl==-1)
          
          class.svm <- ifelse(pi.svm.result[[i]]<=0.5, -1, 1)
          idx.svm <- which(class.svm==-1)
          
          predict.p.fl[[i]][idx.fl] <- 1-pi.fl.result[[i]][idx.fl]
          predict.p.svm[[i]][idx.svm] <- 1-pi.svm.result[[i]][idx.svm]
        }
        
        # calculate the difference
        diff.svm <- rep(0,n.sim)
        diff.fl <- rep(0,n.sim)
        w.diff.svm <- rep(0,n.sim)
        w.diff.fl <- rep(0,n.sim)
        
        ans.p1 <- ans.p
        for (i in 1:n.sim){
          idx1 <- which(c(ans.p1[[i]]-1)==0)
          ans.p1[[i]][idx1] <- c(1-delta)
          
          idx2 <- which(c(ans.p1[[i]]-0)==0)
          ans.p1[[i]][idx2] <- c(delta)
          
          weight <- sqrt(ans.p1[[i]]*(1-ans.p1[[i]]))
          w.diff.svm[i] <- mean(weight*abs(ans.p[[i]] - predict.p.svm[[i]]))
          w.diff.fl[i] <- mean(weight*abs(ans.p[[i]] - predict.p.fl[[i]]))
          
          diff.svm[i] <- mean(abs(ans.p[[i]] - predict.p.svm[[i]]))
          diff.fl[i] <- mean(abs(ans.p[[i]] - predict.p.fl[[i]]))
        } # End of calculating p.diff
        
        a <- round(svm.cre, digits = 3); b <- round(fl.cre, digits = 3)
        c <- round(svm.acc, digits = 3); d <- round(fl.acc, digits = 3)
        e <- round(mean(diff.svm), digits = 3); f <- round(mean(diff.fl), digits = 3)
        g <- round(mean(w.diff.svm), digits = 3); h <- round(mean(w.diff.fl), digits = 3)
        
        rslt.1[j,] <- cbind(a,b,c,d,e,f,g,h)
        # colnames(result) <- c("svm.cre","fl.cre","svm.acc","fl.acc","diff.svm","diff.fl","w.diff.svm","w.diff.fl")
        # rownames(result) <- c("0.5", "1")
      } # End of error loop
      
      result[(r-1)*length(n.train)+trn,] <- rslt.1[j,]
    } # End of training (100~500)
  }
  
  result.by.gen[[n.fn]] <- result
}
